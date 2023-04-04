use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use bstr::ByteSlice;
use chrono::Utc;
use clap::Parser;
use log::{debug, info, warn};
use noodles::core::Position;
use noodles::fasta::fai;
use noodles::fasta::record::{Definition, Sequence};
use noodles::fasta::repository::adapters::IndexedReader;
use noodles::fasta::repository::Adapter;
use noodles::gff::record::Strand;
use noodles::{fasta, gff};
use rayon::prelude::*;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use thiserror::Error;

use drprg::{
    find_prg_index_in, index_fasta, index_vcf, revcomp, Bcftools, GffExt, MakePrg,
    Pandora,
};
use drprg::{MultipleSeqAligner, PathExt};

use crate::cli::check_path_exists;
use crate::config::Config;
use crate::expert::{ExpertRules, RuleExt};
use crate::panel::{Panel, PanelError, PanelExt, PanelRecord};
use crate::Runner;

static META: &str = "##";
const VERSION: &str = env!("CARGO_PKG_VERSION");
static DEFAULT_KMER_SIZE: u32 = 15;
static DEFAULT_WINDOW_SIZE: u32 = 11;
static DEFAULT_PADDING: u32 = 100;
static DEFAULT_MIN_MATCH_LEN: u32 = 5;
static DEFAULT_MAX_NESTING: u32 = 5;

lazy_static! {
    static ref CURRENT_DATE: String = format!("{}", Utc::now().format("%Y%m%d"));
}

#[derive(Parser, Debug, Default)]
pub struct Build {
    /// Path to pandora executable. Will try in src/ext or $PATH if not given
    #[clap(
        short = 'p',
        long = "pandora",
        value_parser,
        hide_short_help = true,
        value_name = "FILE"
    )]
    pandora_exec: Option<PathBuf>,
    /// Path to make_prg executable. Will try in src/ext or $PATH if not given
    #[clap(
        short = 'm',
        long = "makeprg",
        value_parser,
        hide_short_help = true,
        value_name = "FILE"
    )]
    makeprg_exec: Option<PathBuf>,
    /// Path to MAFFT executable. Will try in src/ext or $PATH if not given
    #[clap(
        short = 'M',
        long = "mafft",
        value_parser,
        hide_short_help = true,
        value_name = "FILE"
    )]
    mafft_exec: Option<PathBuf>,
    /// Path to bcftools executable. Will try in src/ext or $PATH if not given
    #[clap(
        short = 'B',
        long = "bcftools",
        value_parser,
        hide_short_help = true,
        value_name = "FILE"
    )]
    bcftools_exec: Option<PathBuf>,
    /// Annotation file that will be used to gather information about genes in catalogue
    #[clap(short = 'a', long = "gff", value_parser = check_path_exists, value_name = "FILE")]
    gff_file: PathBuf,
    /// Panel/catalogue to build index for
    #[clap(short = 'i', long = "panel", value_parser = check_path_exists, value_name = "FILE")]
    panel_file: PathBuf,
    /// Reference genome in FASTA format (must be indexed with samtools faidx)
    #[clap(short = 'f', long = "fasta", value_parser = check_path_exists, value_name = "FILE")]
    reference_file: PathBuf,
    /// An indexed VCF to build the index PRG from. If not provided, then a prebuilt PRG must be
    /// given. See `--prebuilt-prg`
    #[clap(short = 'b', long = "vcf", value_parser = check_path_exists, value_name = "FILE", required_unless_present("prebuilt_prg"))]
    input_vcf: Option<PathBuf>,
    /// Number of bases of padding to add to start and end of each gene
    #[clap(
        short = 'P',
        long = "padding",
        default_value_t = DEFAULT_PADDING,
        value_name = "INT"
    )]
    padding: u32,
    /// Directory to place output
    #[clap(short, long, default_value = ".", value_parser, value_name = "DIR")]
    outdir: PathBuf,
    /// A prebuilt PRG to use.
    ///
    /// Only build the panel VCF and reference sequences - not the PRG. This directory MUST
    /// contain a PRG file named `dr.prg`, along, with a directory called `msas/` that contains an
    /// MSA fasta file for each gene `<gene>.fa`. There can optionally also be a pandora index file, but
    /// if not, the indexing will be performed by drprg. Note: the PRG is expected to contain the
    /// reference sequence for each gene according to the annotation and reference genome given
    /// (along with padding) and must be in the forward strand orientation.
    #[clap(short = 'd', long, value_parser = check_path_exists, value_name = "DIR")]
    prebuilt_prg: Option<PathBuf>,
    /// "Expert rules" to be applied in addition to the catalogue.
    ///
    /// CSV file with blanket rules that describe resistance (or susceptibility). The columns are
    /// <variant type>,<gene>,<start>,<end>,<drug(s)>. See the docs for a detailed explanation.
    #[clap(short, long, value_parser = check_path_exists, value_name = "FILE")]
    rules: Option<PathBuf>,
    /// Minimum number of consecutive characters which must be identical for a match in make_prg
    #[clap(
        short = 'l',
        long = "match-len",
        default_value_t = DEFAULT_MIN_MATCH_LEN,
        hide_short_help = true,
        value_name = "INT"
    )]
    match_len: u32,
    /// Maximum nesting level when constructing the reference graph with make_prg
    #[clap(
        short = 'N',
        long = "max-nesting",
        default_value_t = DEFAULT_MAX_NESTING,
        hide_short_help = true,
        value_name = "INT"
    )]
    max_nesting: u32,
    /// Kmer size to use for pandora
    #[clap(
        short = 'k',
        long,
        default_value_t = DEFAULT_KMER_SIZE,
        hide_short_help = true,
        value_name = "INT"
    )]
    pandora_k: u32,
    /// Window size to use for pandora
    #[clap(
        short = 'w',
        long,
        default_value_t = DEFAULT_WINDOW_SIZE,
        hide_short_help = true,
        value_name = "INT"
    )]
    pandora_w: u32,
    /// Don't index --fasta if an index doesn't exist
    #[clap(short = 'I', long = "no-fai")]
    dont_index_ref: bool,
    /// Don't index --vcf if an index doesn't exist
    #[clap(short = 'C', long = "no-csi")]
    dont_index_vcf: bool,
    /// Version to use for the index
    #[clap(long, default_value_t = CURRENT_DATE.to_string())]
    version: String,
}

impl Build {
    fn load_panel(&self) -> Result<Panel, anyhow::Error> {
        Panel::from_csv(&self.panel_file)
    }

    /// A utility function that checks if a VCF index exists and optionally indexes it
    pub fn check_vcf_index_exists(
        vcf_path: &Path,
        create_index: bool,
    ) -> Result<(), BuildError> {
        for ext in [".csi", ".tbi"] {
            let p = vcf_path.add_extension(ext.as_ref());
            if p.exists() {
                return Ok(());
            }
        }
        if create_index {
            info!("No index found for {:?}. Creating one...", vcf_path);
            index_vcf(vcf_path).map_err(|src| {
                BuildError::IndexError(format!(
                    "Failed to create an index file for {vcf_path:?} due to {src:?}",
                ))
            })
        } else {
            Err(BuildError::IndexError(format!(
                "No index exists for {vcf_path:?} and creation was not requested",
            )))
        }
    }

    /// A utility function that checks if a VCF index exists and optionally indexes it
    pub fn check_fasta_index_exists(
        fasta_path: &Path,
        create_index: bool,
    ) -> Result<(), BuildError> {
        let p = fasta_path.add_extension(".fai".as_ref());
        if p.exists() {
            Ok(())
        } else if create_index {
            info!("No index found for {:?}. Creating one...", fasta_path);
            match index_fasta(fasta_path) {
                Err(src) => Err(BuildError::IndexError(format!(
                    "Failed to create an index file for {fasta_path:?} due to {src:?}",
                ))),
                Ok(_) => Ok(()),
            }
        } else {
            Err(BuildError::IndexError(format!(
                "No index exists for {fasta_path:?} and creation was not requested",
            )))
        }
    }

    fn create_vcf_header(
        &self,
        annotations: &HashMap<String, gff::Record>,
    ) -> bcf::Header {
        let mut header = bcf::Header::new();
        header.push_record(format!("{META}source=drprgV{VERSION}").as_bytes());
        // add contigs to vcf header
        for (gene, gff_record) in annotations {
            let length: u64 = ((usize::from(gff_record.end()) + 1)
                - usize::from(gff_record.start())
                + (&self.padding * 2) as usize) as u64;
            header.push_record(&vcf_contig_field(gene, length));
        }
        for entry in PanelRecord::vcf_header_entries() {
            header.push_record(entry);
        }
        header
    }

    fn rules_path(&self) -> PathBuf {
        self.outdir.join("rules.csv")
    }
    fn prg_path(&self) -> PathBuf {
        self.outdir.join("dr.prg")
    }
    fn msa_dir(&self) -> PathBuf {
        self.outdir.join("msas")
    }
    fn prg_index_path(&self) -> PathBuf {
        let ext = format!("k{}.w{}.idx", self.pandora_k, self.pandora_w);
        self.prg_path().add_extension(ext.as_str().as_ref())
    }
    fn prg_index_kmer_prgs_path(&self) -> PathBuf {
        self.outdir.join("kmer_prgs")
    }
    fn reference_index_file(&self) -> PathBuf {
        self.reference_file.add_extension(".fai".as_ref())
    }
    fn organise_prebuilt_prg(&self) -> Result<(), BuildError> {
        let prebuilt_dir = match &self.prebuilt_prg {
            Some(d) => d.canonicalize().unwrap(), // we know it exists so safe to unwrap
            _ => return Ok(()),
        };
        let prg = prebuilt_dir.join(self.prg_path().file_name().unwrap());
        if !prg.exists() {
            return Err(BuildError::MissingFile(prg));
        }
        let msa_dir = prebuilt_dir.join(self.msa_dir().file_name().unwrap());
        if !msa_dir.exists() {
            return Err(BuildError::MissingFile(msa_dir));
        }
        let prg_index = find_prg_index_in(&prebuilt_dir).ok_or_else(|| {
            BuildError::MissingFile(prebuilt_dir.join("dr.prg.kX.wY.idx"))
        })?;
        let prg_index_kmer_prgs =
            prebuilt_dir.join(self.prg_index_kmer_prgs_path().file_name().unwrap());
        let index_exists = prg_index.exists() && prg_index_kmer_prgs.exists();

        if self.outdir == prebuilt_dir {
            Ok(())
        } else {
            let copyopts = fs_extra::dir::CopyOptions {
                overwrite: true,
                skip_exist: false,
                buffer_size: 64000,
                copy_inside: true,
                content_only: false,
                depth: 0,
            };
            let mut to_copy = vec![prg, msa_dir];
            if index_exists {
                to_copy.push(prg_index);
                to_copy.push(prg_index_kmer_prgs);
            }
            match fs_extra::copy_items(&to_copy, &self.outdir, &copyopts) {
                Ok(_) => Ok(()),
                Err(e) => Err(BuildError::CopyError(format!(
                    "Failed to copy prebuilt PRG files with error message {e}"
                ))),
            }
        }
    }
}

impl Runner for Build {
    fn run(&mut self) -> Result<()> {
        if !self.outdir.exists() {
            info!("Outdir doesn't exist...creating...");
            std::fs::create_dir(&self.outdir)
                .context(format!("Failed to create {:?}", &self.outdir))?;
        }
        self.outdir = self
            .outdir
            .canonicalize()
            .context("Failed to canonicalize outdir")?;

        if self.prebuilt_prg.is_some() {
            self.organise_prebuilt_prg()
                .context("Failed to organise prebuilt PRG files")?;
        } else {
            Build::check_vcf_index_exists(
                self.input_vcf.as_ref().unwrap().as_path(),
                !self.dont_index_vcf,
            )?;
        }

        Build::check_fasta_index_exists(&self.reference_file, !self.dont_index_ref)?;

        let premsa_dir = self.outdir.join("premsa");
        let msa_dir = self.msa_dir();
        info!("Building panel index...");

        info!("Loading the panel...");
        let panel: Panel = self.load_panel()?;
        let mut genes: HashSet<String> = panel.keys().map(|k| k.to_owned()).collect();
        info!("Loaded panel");

        if let Some(p) = &self.rules {
            let rules =
                ExpertRules::from_csv(p).context("Failed to read expert rules")?;
            for g in rules.keys() {
                let _ = genes.insert(g.to_owned());
            }
            fs::copy(p, self.rules_path())
                .context("Failed to copy expert rules file")?;
        }

        info!("Loading genome annotation for panel genes...");
        debug!("Panel genes: {:?}", genes);
        let gff_reader = File::open(&self.gff_file)
            .map(BufReader::new)
            .map(gff::Reader::new)?;
        let annotations = load_annotations_for_genes(gff_reader, &genes);
        if genes.len() != annotations.len() {
            let annotation_genes: HashSet<String> =
                annotations.keys().cloned().collect();
            let diff = genes.difference(&annotation_genes);
            warn!("Couldn't load annotations for genes {:?}", diff);
        }
        info!("Loaded annotations");

        let panel_vcf_path = self.outdir.join("panel.bcf");
        let gene_refs_path = self.outdir.join("genes.fa");
        debug!("Creating VCF header...");
        let vcf_header = self.create_vcf_header(&annotations);
        debug!("VCF header created");
        let unsorted_panel_vcf_path = panel_vcf_path.with_extension("unsorted.bcf");
        {
            let mut vcf_writer = bcf::Writer::from_path(
                &unsorted_panel_vcf_path,
                &vcf_header,
                false,
                bcf::Format::Bcf,
            )?;
            debug!("Loading the reference genome index...");

            let fa_reader = File::open(&self.reference_file)
                .map(BufReader::new)
                .map(fasta::Reader::new)?;
            let index = File::open(self.reference_index_file())
                .map(BufReader::new)
                .map(fai::Reader::new)?
                .read_index()?;
            let mut faidx = IndexedReader::new(fa_reader, index);
            debug!("Loaded the reference genome index");
            let mut fa_writer =
                File::create(&gene_refs_path).map(BufWriter::new).map(|f| {
                    fasta::Writer::builder(f)
                        .set_line_base_count(usize::MAX)
                        .build()
                })?;
            if !premsa_dir.exists() {
                debug!("Pre-MSA directory doesn't exist...creating...");
                std::fs::create_dir(&premsa_dir)
                    .context(format!("Failed to create {:?}", &premsa_dir))?;
            } else {
                debug!("Existing pre-MSA directory found...removing...");
                std::fs::remove_dir_all(&premsa_dir)
                    .and_then(|_| std::fs::create_dir(&premsa_dir))
                    .context(format!(
                        "Failed to remove and recreate {:?}",
                        &premsa_dir
                    ))?;
            }

            info!("Converting the panel to a VCF...");
            for (gene, gff_record) in &annotations {
                let seq =
                    extract_gene_from_index(gff_record, &mut faidx, self.padding)?;
                let definition =
                    Definition::new(gene, Some(format!("padding={}", self.padding)));
                let fa_record =
                    fasta::Record::new(definition, Sequence::from(seq.to_owned()));
                fa_writer.write_record(&fa_record)?;

                let premsa_path = premsa_dir.join(format!("{gene}.fa"));
                let mut premsa_writer =
                    File::create(premsa_path).map(BufWriter::new).map(|f| {
                        fasta::Writer::builder(f)
                            .set_line_base_count(usize::MAX)
                            .build()
                    })?;
                premsa_writer.write_record(&fa_record)?;

                let panel_records_for_gene = &panel[gene];
                for panel_record in panel_records_for_gene {
                    let mut vcf_record = vcf_writer.empty_record();
                    match panel_record.to_vcf(&mut vcf_record, &seq, self.padding) {
                        Err(PanelError::PosOutOfRange(pos, gene)) => {
                            warn!(
                                "Position {} is out of range for gene {} [Skipping]",
                                pos, gene
                            );
                            continue;
                        }
                        Err(e) => return Err(anyhow!(e)),
                        _ => (),
                    }

                    let strand = gff_record.strand();
                    vcf_record
                        .push_info_string(b"ST", &[strand.as_ref().as_bytes()])
                        .context(format!(
                            "Couldn't set INFO field ST for {}",
                            panel_record.name()
                        ))?;
                    vcf_writer.write(&vcf_record)?;
                }
            }
        }
        debug!("Indexing the genes reference fasta...");
        index_fasta(&gene_refs_path)
            .context("Failed to index genes reference fasta")?;

        debug!("Sorting the original panel VCF...");
        let bcftools = Bcftools::from_path(&self.bcftools_exec)?;
        bcftools
            .sort(&unsorted_panel_vcf_path, &panel_vcf_path)
            .context("Failed to sort the panel VCF")?;
        if let Err(e) = std::fs::remove_file(unsorted_panel_vcf_path) {
            warn!(
                "Failed to remove the unsorted/original panel VCF with error:\n{}",
                e.to_string()
            );
        }
        debug!("Panel VCF sorted and original removed");
        info!(
            "Panel successfully converted to a VCF at {:?}",
            panel_vcf_path
        );

        info!("Indexing the panel VCF...");
        index_vcf(&panel_vcf_path).context("Failed to index the panel VCF")?;
        info!("Panel VCF successfully indexed");

        if self.prebuilt_prg.is_none() {
            info!("Applying input VCF variants to reference (pre-MSA)");
            // we can unwrap as it MUST be some if prebuilt prg is none
            let input_vcf_path = self.input_vcf.as_ref().unwrap();
            let input_vcf = bcf::Reader::from_path(input_vcf_path)
                .context("Failed to open input VCF file")?;
            let samples = input_vcf.header().samples();

            let consensus_dir = self.outdir.join("consensus");
            if !consensus_dir.exists() {
                debug!("Consensus directory doesn't exist...creating...");
                std::fs::create_dir(&consensus_dir)
                    .context(format!("Failed to create {:?}", &consensus_dir))?;
            } else {
                debug!("Existing consensus directory found...removing...");
                std::fs::remove_dir_all(&consensus_dir)
                    .and_then(|_| std::fs::create_dir(&consensus_dir))
                    .context(format!(
                        "Failed to remove and recreate {:?}",
                        &consensus_dir
                    ))?;
            }

            debug!("Writing consensus sequence reference file...");
            let consensus_ref_path = consensus_dir.join("ref.fa");
            let mut fa_reader = File::open(&gene_refs_path)
                .map(BufReader::new)
                .map(fasta::Reader::new)?;
            {
                let mut fa_writer = File::create(&consensus_ref_path)
                    .map(BufWriter::new)
                    .map(|f| {
                        fasta::Writer::builder(f)
                            .set_line_base_count(usize::MAX)
                            .build()
                    })?;
                for result in fa_reader.records() {
                    let record = result?;
                    let gene = record.definition().name();
                    // safe to unwrap as we built this fasta from the annotations
                    let gff_record = annotations.get(gene).unwrap();
                    let seq = match gff_record.strand() {
                        Strand::Reverse => Ok(revcomp(record.sequence().as_ref())),
                        Strand::Forward => Ok(record.sequence().as_ref().to_vec()),
                        _ => Err(BuildError::MissingStrand(gene.to_string())),
                    }?;
                    let out_record = fasta::Record::new(
                        Definition::new(gene, None),
                        Sequence::from(seq),
                    );
                    fa_writer
                        .write_record(&out_record)
                        .context("Failed to write consensus reference sequence")?;
                }
            }

            samples.into_par_iter().try_for_each(|sample| {
                let s = sample.to_str_lossy();
                let consensus_path = consensus_dir.join(format!("{s}.fa"));
                let args = &[
                    "-H",
                    "A",
                    "-s",
                    &s,
                    "-f",
                    &consensus_ref_path.to_string_lossy(),
                ];
                bcftools.consensus(input_vcf_path, &consensus_path, args)
            })?;
            debug!(
                "Generated consensus sequences for all {} sample(s) in the input VCF",
                input_vcf.header().sample_count()
            );
            std::fs::remove_file(&consensus_ref_path)
                .context("Failed to remove consensus reference file")?;

            {
                debug!("Combining consensus sequences into pre-MSAs...");
                for sample in input_vcf.header().samples() {
                    let s = sample.to_str_lossy().to_string();
                    let sample_consensus_path = consensus_dir.join(format!("{s}.fa"));
                    let mut sample_consensus_fasta = File::open(sample_consensus_path)
                        .map(BufReader::new)
                        .map(fasta::Reader::new)?;
                    for result in sample_consensus_fasta.records() {
                        let consensus_record = result.context(format!("Couldn't parse fasta record from bcftools consensus for sample {}", &s))?;
                        let gene = consensus_record.name();
                        let is_reverse =
                            annotations.get(gene).unwrap().strand() == Strand::Reverse;
                        let gene_pre_msa_path = premsa_dir.join(format!("{gene}.fa"));
                        {
                            let mut premsa_writer = File::options()
                                .append(true)
                                .open(gene_pre_msa_path)
                                .map(BufWriter::new)
                                .map(|f| {
                                    fasta::Writer::builder(f)
                                        .set_line_base_count(usize::MAX)
                                        .build()
                                })?;
                            let def = Definition::new(&s, None);
                            let seq = if is_reverse {
                                revcomp(consensus_record.sequence().as_ref())
                            } else {
                                consensus_record.sequence().as_ref().to_vec()
                            };
                            let out_record =
                                fasta::Record::new(def, Sequence::from(seq));
                            premsa_writer.write_record(&out_record)?;
                        }
                    }
                }
            }

            info!("Generating multiple sequence alignments and reference graphs for all genes and their variants...");
            let mafft = MultipleSeqAligner::from_path(&self.mafft_exec)?;

            if !msa_dir.exists() {
                debug!("MSA directory doesn't exist...creating...");
                std::fs::create_dir(&msa_dir)
                    .context(format!("Failed to create {:?}", &msa_dir))?;
            } else {
                info!("Existing MSA directory found...removing...");
                std::fs::remove_dir_all(&msa_dir)
                    .and_then(|_| std::fs::create_dir(&msa_dir))
                    .context(format!("Failed to remove and recreate {:?}", &msa_dir))?;
            }

            genes.par_iter().try_for_each(|gene| {
                let premsa_path = premsa_dir.join(format!("{gene}.fa"));
                let msa_path = msa_dir.join(format!("{gene}.fa"));
                debug!("Running MSA for {}", gene);
                mafft.run_with(
                    &premsa_path,
                    &msa_path,
                    ["--auto", "--thread", "-1"],
                    true,
                )
            })?;
            info!("Successfully generated MSAs");

            info!("Building reference graphs for genes...");
            let makeprg = MakePrg::from_path(&self.makeprg_exec)?;
            let make_prg_args = &[
                "-t",
                &rayon::current_num_threads().to_string(),
                "-L",
                &self.match_len.to_string(),
                "-N",
                &self.max_nesting.to_string(),
            ];
            makeprg.from_msas_with(&msa_dir, &self.prg_path(), make_prg_args)?;
            info!("Successfully created panel reference graph");
        }
        info!("Indexing reference graph with pandora...");
        let pandora = Pandora::from_path(&self.pandora_exec)?;
        let pandora_args = &[
            "-t",
            &rayon::current_num_threads().to_string(),
            "-k",
            &self.pandora_k.to_string(),
            "-w",
            &self.pandora_w.to_string(),
        ];
        if self.prebuilt_prg.is_some() && self.prg_index_path().exists() {
            info!("Existing pandora index found...skipping...");
        } else {
            pandora.index_with(&self.prg_path(), pandora_args)?;
            info!("Reference graph indexed successfully");
        }
        info!("Panel index built");

        let config_file = &self.outdir.join(".config.toml");
        {
            let mut f = File::create(config_file).map(BufWriter::new)?;
            let config: Config = Config {
                min_match_len: self.match_len,
                max_nesting: self.max_nesting,
                k: self.pandora_k,
                w: self.pandora_w,
                padding: self.padding,
                version: self.version.to_owned(),
            };
            let toml = toml::to_string(&config)?;
            writeln!(f, "{toml}")?;
        }

        debug!("Cleaning up temporary files...");
        if premsa_dir.exists() {
            std::fs::remove_dir_all(&premsa_dir)
                .context(format!("Failed to remove {:?}", &msa_dir))?;
        }
        Ok(())
    }
}

/// A collection of custom errors relating to the build component of this package
#[derive(Error, Debug, PartialEq, Eq)]
pub enum BuildError {
    /// Contig is missing from the FASTA index
    #[error("Contig {0} is missing from the reference genome index")]
    MissingContig(String),
    /// Couldn't fetch a gene
    #[error("Couldn't fetch contig {0} with start {1} and end {2}")]
    FetchError(String, u64, u64),
    /// GFF record doesn't have a strand
    #[error("No strand information in GFF for {0}")]
    MissingStrand(String),
    /// An expected file is missing
    #[error("Missing expected file {0:?}")]
    MissingFile(PathBuf),
    /// File copying error
    #[error("{0}")]
    CopyError(String),
    /// An issue with an index file
    #[error("{0}")]
    IndexError(String),
}

fn load_annotations_for_genes<R>(
    mut reader: gff::Reader<R>,
    genes: &HashSet<String>,
) -> HashMap<String, gff::Record>
where
    R: BufRead,
{
    let mut records = reader.records();
    let mut annotations: HashMap<String, gff::Record> = HashMap::new();

    while let Some(Ok(record)) = records.next() {
        if record.ty() != "gene" {
            continue;
        }

        if let Some(entry) = record
            .attributes()
            .iter()
            .find(|&entry| entry.key() == "Name")
        {
            if genes.contains(entry.value()) {
                annotations.insert(entry.value().to_string(), record.to_owned());
            }
        }
    }
    annotations
}

fn extract_gene_from_index<R>(
    record: &gff::Record,
    faidx: &mut IndexedReader<R>,
    padding: u32,
) -> Result<Vec<u8>, BuildError>
where
    R: BufRead + Seek,
{
    let idx_record = faidx.get(record.reference_sequence_name());
    let contig = match &idx_record {
        Some(Ok(contig)) => contig,
        _ => {
            return Err(BuildError::MissingContig(
                record.reference_sequence_name().to_string(),
            ))
        }
    };
    let contig_len = contig.sequence().len() as u64;

    // GFF3 start is 1-based, inclusive. We want to make it 0-based, inclusive
    let start: u64 = max(
        usize::from(record.start()) as isize - padding as isize - 1,
        0,
    ) as u64;
    // GFF3 end is 1-based, inclusive. We want to make it 0-based, exclusive
    let end: u64 = min(
        usize::from(record.end()) as u64 + padding as u64,
        contig_len,
    );
    if start > end {
        return Err(BuildError::FetchError(
            record.reference_sequence_name().to_string(),
            start,
            end,
        ));
    }

    // noodles sequence slice/positions are 1-based
    let fetch_slice = Position::try_from((start + 1) as usize).unwrap()
        ..=Position::try_from(end as usize).unwrap();
    match contig.sequence().slice(fetch_slice) {
        None => Err(BuildError::FetchError(
            record.reference_sequence_name().to_string(),
            start,
            end,
        )),
        Some(seq) => {
            let is_rev = match record.strand() {
                Strand::Reverse => true,
                Strand::Forward => false,
                _ => {
                    return Err(BuildError::MissingStrand(
                        record.name().unwrap_or("No Name").to_string(),
                    ))
                }
            };
            if is_rev {
                Ok(revcomp(seq.as_ref()))
            } else {
                Ok(Vec::from(seq.as_ref()))
            }
        }
    }
}

fn vcf_contig_field(id: &str, length: u64) -> Vec<u8> {
    format!("{META}contig=<ID={id},length={length}>")
        .as_bytes()
        .to_owned()
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Read, Write};
    use std::iter::FromIterator;
    use std::str::FromStr;

    use crate::panel::{Residue, Variant};

    use super::*;

    const PANEL: &str = "tests/cases/panel.tsv";
    const ANNOTATION: &str = "tests/cases/ann.gff3";
    const VCF: &str = "tests/cases/build/input.bcf";
    const REF: &str = "tests/cases/ref.fa";

    #[test]
    fn load_annotations_when_no_genes_in_common_returns_empty() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let reader = gff::Reader::new(GFF);
        let genes = HashSet::from_iter(vec!["geneX".to_string()]);

        let actual = load_annotations_for_genes(reader, &genes);
        assert!(actual.is_empty())
    }

    #[test]
    fn load_annotations_for_genes_one_gene_in_common() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let reader = gff::Reader::new(GFF);
        let genes = HashSet::from_iter(vec!["geneX".to_string(), "dnaA".to_string()]);

        let actual = load_annotations_for_genes(reader, &genes);
        let (actual_gene, actual_record) = actual.iter().next().unwrap();
        assert_eq!(actual_gene, "dnaA");
        assert_eq!(actual_record.end(), Position::new(1524).unwrap())
    }

    #[test]
    fn load_annotations_for_genes_is_cds_returns_empty() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tCDS\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let reader = gff::Reader::new(GFF);
        let genes = HashSet::from_iter(vec!["geneX".to_string(), "dnaA".to_string()]);

        let actual = load_annotations_for_genes(reader, &genes);
        assert!(actual.is_empty())
    }

    #[test]
    fn extract_gene_from_index_gene_not_in_index() {
        let padding = 0;
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tCDS\t1\t1524\t.\t+\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";

        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap_err();
        let expected = BuildError::MissingContig("NC_000962.3".to_string());
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_interval_out_of_bounds() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t100\t1524\t.\t+\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap_err();
        let expected = BuildError::FetchError("chr1".to_string(), 99, 16);
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_first_base() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t1\t.\t+\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"G";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_too_much_padding_left_wraps_to_start() {
        let padding = 2;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t1\t.\t+\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"GTA";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_too_much_padding_right_wraps_to_end() {
        let padding = 4;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t16\t16\t.\t+\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"ACCCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_padding_start_and_end_exactly_the_same_as_gene() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t16\t.\t+\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"GTAGGCTGAAAACCCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_on_reverse_strand() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t16\t.\t-\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = revcomp(b"GTAGGCTGAAAACCCC");
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_strand() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t16\t.\t.\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap_err();
        let expected = BuildError::MissingStrand("dnaA".to_string());
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_padding_start_same_as_gene_end_minus_one_from_gene_length(
    ) {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t15\t.\t+\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"GTAGGCTGAAAACCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_padding_start_plus_one_from_gene_start_end_same_as_gene(
    ) {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t2\t16\t.\t+\t0\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let fa_reader = fasta::Reader::new(Cursor::new(FASTA_FILE));
        let index = fai::Reader::new(FAI_FILE).read_index().unwrap();
        let mut faidx = IndexedReader::new(fa_reader, index);

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"TAGGCTGAAAACCCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn test_vcf_contig_field() {
        let id = "chr1";
        let length = 33;

        let actual = vcf_contig_field(id, length);
        let expected = b"##contig=<ID=chr1,length=33>";

        assert_eq!(actual, expected)
    }

    #[test]
    fn load_panel_single_record() {
        let gene = "pncA";
        let variant = "G6T";
        let residue = "DNA";
        let drugs = "Drug1";
        let contents: &str = &[gene, variant, residue, drugs].join("\t");
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        let builder = Build {
            panel_file: PathBuf::from(file.path()),
            ..Default::default()
        };

        let actual = builder.load_panel().unwrap();
        let expected = HashMap::from_iter(vec![(
            gene.to_string(),
            HashSet::from_iter(vec![PanelRecord {
                gene: gene.to_string(),
                variant: Variant::from_str(variant).unwrap(),
                residue: Residue::from_str(residue).unwrap(),
                drugs: HashSet::from_iter(vec![drugs.to_string()]),
            }]),
        )]);

        assert_eq!(actual, expected)
    }

    #[test]
    fn load_panel_duplicate_record() {
        let gene = "pncA";
        let variant = "G6T";
        let residue = "DNA";
        let drugs = "Drug1";
        let contents: &str = &[gene, variant, residue, drugs].join("\t");
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        let builder = Build {
            panel_file: PathBuf::from(file.path()),
            ..Default::default()
        };

        let actual = builder.load_panel().unwrap();
        let expected = HashMap::from_iter(vec![(
            gene.to_string(),
            HashSet::from_iter(vec![PanelRecord {
                gene: gene.to_string(),
                variant: Variant::from_str(variant).unwrap(),
                residue: Residue::from_str(residue).unwrap(),
                drugs: HashSet::from_iter(vec![drugs.to_string()]),
            }]),
        )]);

        assert_eq!(actual, expected)
    }

    #[test]
    fn load_panel_wrong_delimiter() {
        let gene = "pncA";
        let variant = "G6T";
        let residue = "DNA";
        let drugs = "Drug1";
        let contents: &str = &[gene, variant, residue, drugs].join(",");
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        let builder = Build {
            panel_file: PathBuf::from(file.path()),
            ..Default::default()
        };

        let actual = builder.load_panel();
        assert!(actual.is_err())
    }

    #[test]
    fn load_panel_has_header() {
        let header: &str = &["gene", "variant", "residue", "drugs"].join("\t");
        let gene = "pncA";
        let variant = "G6T";
        let residue = "DNA";
        let drugs = "Drug1";
        let contents: &str = &[gene, variant, residue, drugs].join("\t");
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(header.as_bytes()).unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        let builder = Build {
            panel_file: PathBuf::from(file.path()),
            ..Default::default()
        };

        let actual = builder.load_panel();
        assert!(actual.is_err())
    }

    #[test]
    fn load_panel_path_doesnt_exist() {
        let builder = Build {
            panel_file: PathBuf::from("foobar"),
            ..Default::default()
        };

        let actual = builder.load_panel();
        assert!(actual.is_err())
    }

    #[test]
    fn create_vcf_header() {
        let builder = Build::default();
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Name=dnaA;gene=dnaA\n";
        let reader = gff::Reader::new(GFF);
        let genes = HashSet::from_iter(vec!["dnaA".to_string()]);
        let annotations = load_annotations_for_genes(reader, &genes);
        let tmp = tempfile::NamedTempFile::new().unwrap();

        let header = builder.create_vcf_header(&annotations);

        let writer =
            bcf::Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf)
                .unwrap();
        let view = writer.header();

        assert_eq!(view.contig_count(), 1);
        assert_eq!(view.rid2name(0).unwrap(), b"dnaA")
    }

    #[test]
    fn build_runner() {
        let outdir = tempfile::tempdir().unwrap();
        let mut builder = Build {
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            makeprg_exec: Some(PathBuf::from("src/ext/make_prg")),
            mafft_exec: Some(PathBuf::from("src/ext/mafft/bin/mafft")),
            bcftools_exec: Some(PathBuf::from("src/ext/bcftools")),
            gff_file: ANNOTATION.parse().unwrap(),
            panel_file: PANEL.parse().unwrap(),
            input_vcf: Some(VCF.parse().unwrap()),
            reference_file: REF.parse().unwrap(),
            padding: 100,
            outdir: PathBuf::from(outdir.path()),
            match_len: 5,
            max_nesting: 7,
            pandora_w: 14,
            pandora_k: 15,
            ..Build::default()
        };
        let result = builder.run();
        assert!(result.is_ok());

        let mut file1 = File::open("tests/cases/expected/dr.prg").unwrap();
        let mut file2 =
            File::open(format!("{}/dr.prg", outdir.path().to_string_lossy())).unwrap();

        let mut contents = String::new();
        file1.read_to_string(&mut contents).unwrap();
        let mut other = String::new();
        file2.read_to_string(&mut other).unwrap();

        let mut sorted1 = contents.as_bytes().to_owned();
        sorted1.sort_unstable();

        let mut sorted2 = other.as_bytes().to_owned();
        sorted2.sort_unstable();

        assert_eq!(sorted1, sorted2);
    }

    #[test]
    fn build_runner_outdir_doesnt_exist() {
        let outdir = std::env::temp_dir().join("foooo");
        if outdir.exists() {
            std::fs::remove_dir_all(&outdir).unwrap();
        }
        let mut builder = Build {
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            makeprg_exec: Some(PathBuf::from("src/ext/make_prg")),
            mafft_exec: Some(PathBuf::from("src/ext/mafft/bin/mafft")),
            bcftools_exec: Some(PathBuf::from("src/ext/bcftools")),
            gff_file: ANNOTATION.parse().unwrap(),
            panel_file: PANEL.parse().unwrap(),
            input_vcf: Some(VCF.parse().unwrap()),
            reference_file: REF.parse().unwrap(),
            padding: 100,
            outdir: outdir.to_owned(),
            match_len: 5,
            max_nesting: 7,
            pandora_w: 14,
            pandora_k: 15,
            ..Build::default()
        };
        let result = builder.run();
        assert!(result.is_ok());

        let mut file1 = File::open("tests/cases/expected/dr.prg").unwrap();
        let mut file2 =
            File::open(format!("{}/dr.prg", outdir.to_string_lossy())).unwrap();

        let mut contents = String::new();
        file1.read_to_string(&mut contents).unwrap();
        let mut other = String::new();
        file2.read_to_string(&mut other).unwrap();

        let mut sorted1 = contents.as_bytes().to_owned();
        sorted1.sort_unstable();

        let mut sorted2 = other.as_bytes().to_owned();
        sorted2.sort_unstable();

        assert_eq!(sorted1, sorted2);
    }

    #[test]
    fn check_vcf_and_index_exists_it_doesnt() {
        let p = PathBuf::from("tests/cases/predict/in.vcf");
        let actual = Build::check_vcf_index_exists(&p, false);
        assert!(actual.is_err())
    }

    #[test]
    fn check_vcf_and_index_exists_it_does() {
        let path = PathBuf::from("tests/cases/build/input.bcf");
        let actual = Build::check_vcf_index_exists(&path, false);
        assert!(actual.is_ok())
    }
}
