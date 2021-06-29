use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::io::{Read, Seek};
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};
use bio::alphabets::dna::revcomp;
use bio::io::{fasta, gff};
use log::{debug, info, warn};
use rayon::prelude::*;
use rust_htslib::bcf;
use structopt::clap::AppSettings;
use structopt::StructOpt;
use thiserror::Error;

use drprg::{index_vcf, Bcftools, MakePrg, Pandora};
use drprg::{MultipleSeqAligner, PathExt};

use crate::cli::check_path_exists;
use crate::panel::{Panel, PanelError, PanelExt, PanelRecord};
use crate::Runner;

static META: &str = "##";
const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(StructOpt, Debug, Default)]
#[structopt(setting = AppSettings::DeriveDisplayOrder)]
pub struct Build {
    /// Path to pandora executable. Will try in src/ext or $PATH if not given
    #[structopt(
        short = "p",
        long = "pandora",
        parse(from_os_str),
        hidden_short_help = true,
        value_name = "FILE"
    )]
    pandora_exec: Option<PathBuf>,
    /// Path to make_prg executable. Will try in src/ext or $PATH if not given
    #[structopt(
        short = "m",
        long = "makeprg",
        parse(from_os_str),
        hidden_short_help = true,
        value_name = "FILE"
    )]
    makeprg_exec: Option<PathBuf>,
    /// Path to MAFFT executable. Will try in src/ext or $PATH if not given
    #[structopt(
        short = "M",
        long = "mafft",
        parse(from_os_str),
        hidden_short_help = true,
        value_name = "FILE"
    )]
    mafft_exec: Option<PathBuf>,
    /// Path to bcftools executable. Will try in src/ext or $PATH if not given
    #[structopt(
        short = "B",
        long = "bcftools",
        parse(from_os_str),
        hidden_short_help = true,
        value_name = "FILE"
    )]
    bcftools_exec: Option<PathBuf>,
    /// Annotation file that will be used to gather information about genes in panel
    #[structopt(short = "a", long = "gff", parse(try_from_os_str = check_path_exists), value_name = "FILE")]
    gff_file: PathBuf,
    /// Panel to build index from
    #[structopt(short = "i", long = "panel", parse(try_from_os_str = check_path_exists), value_name = "FILE")]
    panel_file: PathBuf,
    /// Reference genome in FASTA format (must be indexed with samtools faidx)
    #[structopt(short = "f", long = "fasta", parse(try_from_os_str = check_path_exists), value_name = "FILE")]
    reference_file: PathBuf,
    /// Number of bases of padding to add to start and end of each gene
    #[structopt(
        short = "P",
        long = "padding",
        default_value = "100",
        value_name = "INT"
    )]
    padding: u32,
    /// Directory to place output
    #[structopt(
        short,
        long,
        default_value = ".",
        parse(from_os_str),
        value_name = "DIR"
    )]
    outdir: PathBuf,
    /// Minimum number of consecutive characters which must be identical for a match in make_prg
    #[structopt(
        short = "-l",
        long = "--match-len",
        default_value = "7",
        hidden_short_help = true,
        value_name = "INT"
    )]
    match_len: u32,
    /// Maximum nesting level when constructing the panel graph with make_prg
    #[structopt(
        short = "-N",
        long = "--max-nesting",
        default_value = "5",
        hidden_short_help = true,
        value_name = "INT"
    )]
    max_nesting: u32,
    /// Force overwriting existing files. Use this if you want to build from scratch
    #[structopt(short = "F", long)]
    force: bool,
}

impl Build {
    fn load_panel(&self) -> Result<Panel, anyhow::Error> {
        Panel::from_csv(&self.panel_file)
    }

    fn create_vcf_header(
        &self,
        annotations: &HashMap<String, gff::Record>,
    ) -> bcf::Header {
        let mut header = bcf::Header::new();
        header.push_record(format!("{}source=drprgV{}", META, VERSION).as_bytes());
        // add contigs to vcf header
        for (gene, gff_record) in annotations {
            let length: u64 = (gff_record.end() + 1) - gff_record.start()
                + (&self.padding * 2) as u64;
            header.push_record(&*vcf_contig_field(gene, length));
        }
        for entry in PanelRecord::vcf_header_entries() {
            header.push_record(entry);
        }
        header
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
        let premsa_dir = self.outdir.join("premsa");
        let msa_dir = self.outdir.join("msa");
        info!("Building panel index...");

        info!("Loading the panel...");
        let panel: Panel = self.load_panel()?;
        let genes: HashSet<String> = panel.keys().map(|k| k.to_owned()).collect();
        info!("Loaded panel");

        info!("Loading genome annotation for panel genes...");
        debug!("Panel genes: {:?}", genes);
        let gff_reader = gff::Reader::from_file(&self.gff_file, gff::GffType::GFF3)?;
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
        if !self.force && panel_vcf_path.exists() && gene_refs_path.exists() {
            info!("Using existing panel VCF and gene reference files")
        } else {
            debug!("Creating VCF header...");
            let vcf_header = self.create_vcf_header(&annotations);
            debug!("VCF header created");
            let unsorted_panel_vcf_path = panel_vcf_path.with_extension("unsorted.bcf");
            {
                let mut vcf_writer = bcf::Writer::from_path(
                    &unsorted_panel_vcf_path,
                    &vcf_header,
                    false,
                    bcf::Format::BCF,
                )?;
                debug!("Loading the reference genome index...");
                let mut faidx = fasta::IndexedReader::from_file(&self.reference_file)?;
                debug!("Loaded the reference genome index");
                let mut fa_writer = fasta::Writer::to_file(gene_refs_path)?;
                if !premsa_dir.exists() {
                    debug!("Pre-MSA directory doesn't exist...creating...");
                    std::fs::create_dir(&premsa_dir)
                        .context(format!("Failed to create {:?}", &premsa_dir))?;
                } else if self.force {
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
                        extract_gene_from_index(&gff_record, &mut faidx, self.padding)?;
                    fa_writer.write(
                        gene,
                        Some(&format!("padding={}", self.padding)),
                        &seq,
                    )?;
                    let premsa_path = premsa_dir.join(format!("{}.fa", gene));
                    let mut premsa_writer = fasta::Writer::to_file(premsa_path)?;
                    premsa_writer.write(
                        gene,
                        Some(&format!("padding={}", self.padding)),
                        &seq,
                    )?;

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
                        // we use unwrap as we have already confirmed above that the strand is + or -
                        let strand =
                            gff_record.strand().unwrap().strand_symbol().to_owned();
                        vcf_record
                            .push_info_string(b"ST", &[strand.as_bytes()])
                            .context(format!(
                                "Couldn't set INFO field ST for {}",
                                panel_record.name()
                            ))?;
                        vcf_writer.write(&vcf_record)?;

                        let s: usize = vcf_record.pos() as usize;
                        let e: usize = s + vcf_record.inner().rlen as usize;
                        for (i, allele) in vcf_record.alleles()[1..].iter().enumerate()
                        {
                            let mutated_seq =
                                [seq[..s].to_vec(), allele.to_vec(), seq[e..].to_vec()]
                                    .concat();
                            premsa_writer.write(
                                &format!("{}_{}", panel_record.name(), i),
                                None,
                                &mutated_seq,
                            )?;
                        }
                    }
                }
            }
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
        }

        let panel_vcf_index_path = panel_vcf_path.add_extension(".csi".as_ref());
        if !self.force && panel_vcf_index_path.exists() {
            info!("Using existing panel VCF index")
        } else {
            info!("Indexing the panel VCF...");
            index_vcf(&panel_vcf_path).context("Failed to index the panel VCF")?;
            info!("Panel VCF successfully indexed");
        }

        info!("Generating multiple sequence alignments and reference graphs for all genes and their variants...");
        let mafft = MultipleSeqAligner::from_path(&self.mafft_exec)?;

        if !msa_dir.exists() {
            debug!("MSA directory doesn't exist...creating...");
            std::fs::create_dir(&msa_dir)
                .context(format!("Failed to create {:?}", &msa_dir))?;
        } else if self.force {
            info!("Existing MSA directory found...removing...");
            std::fs::remove_dir_all(&msa_dir)
                .and_then(|_| std::fs::create_dir(&msa_dir))
                .context(format!("Failed to remove and recreate {:?}", &msa_dir))?;
        }

        genes.par_iter().try_for_each(|gene| {
            let premsa_path = premsa_dir.join(format!("{}.fa", gene));
            let msa_path = msa_dir.join(format!("{}.fa", gene));
            if !self.force && msa_path.exists() {
                Ok(())
            } else {
                debug!("Running MSA for {}", gene);
                mafft.run_with(&premsa_path, &msa_path, &["--auto", "--thread", "-1"])
            }
        })?;
        info!("Successfully generated MSAs");

        info!("Building reference graphs for genes...");
        let makeprg = MakePrg::from_path(&self.makeprg_exec)?;
        let prg_path = self.outdir.join("dr.prg");
        let prg_update_prgs = self.outdir.join("dr_prgs");
        let prg_update_ds = prg_path.with_extension("update_DS");
        if !self.force
            && prg_path.exists()
            && prg_update_ds.exists()
            && prg_update_prgs.exists()
        {
            info!("Existing reference graph found...skipping...");
        } else {
            let make_prg_args = &[
                "-t",
                &rayon::current_num_threads().to_string(),
                "--min_match_length",
                &self.match_len.to_string(),
                "--max_nesting",
                &self.max_nesting.to_string(),
            ];
            makeprg.from_msas_with(
                &msa_dir,
                &prg_path,
                &prg_update_ds,
                &prg_update_prgs,
                make_prg_args,
            )?;
            info!("Successfully created panel reference graph");
        }
        info!("Indexing reference graph with pandora...");
        let pandora = Pandora::from_path(&self.pandora_exec)?;
        let pandora_args = &["-t", &rayon::current_num_threads().to_string()];
        let pandora_index = prg_path.add_extension(".k15.w14.idx".as_ref());
        if !self.force && pandora_index.exists() {
            info!("Existing pandora index found...skipping...");
        } else {
            pandora.index_with(&prg_path, pandora_args)?;
            info!("Reference graph indexed successfully");
        }
        info!("Panel index built");
        Ok(())
    }
}

/// A collection of custom errors relating to the build component of this package
#[derive(Error, Debug, PartialEq)]
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
}

fn load_annotations_for_genes<R>(
    mut reader: gff::Reader<R>,
    genes: &HashSet<String>,
) -> HashMap<String, gff::Record>
where
    R: std::io::Read,
{
    let mut records = reader.records();
    let mut annotations: HashMap<String, gff::Record> = HashMap::new();

    while let Some(Ok(record)) = records.next() {
        if record.feature_type() != "gene" {
            continue;
        }

        if let Some(g) = record.attributes().get("gene") {
            if genes.contains(g) {
                annotations.insert(g.to_string(), record.to_owned());
            }
        }
    }
    annotations
}

fn extract_gene_from_index<R>(
    record: &gff::Record,
    faidx: &mut fasta::IndexedReader<R>,
    padding: u32,
) -> Result<Vec<u8>, BuildError>
where
    R: Read + Seek,
{
    let mut contig_len: u64 = 0;
    for contig in faidx.index.sequences() {
        if contig.name == record.seqname() {
            contig_len = contig.len;
            break;
        }
    }
    if contig_len == 0 {
        return Err(BuildError::MissingContig(record.seqname().to_string()));
    }

    // GFF3 start is 1-based, inclusive. We want to make it 0-based, inclusive
    let start: u64 = max(
        (*record.start() as isize - padding as isize - 1) as isize,
        0,
    ) as u64;
    // GFF3 end is 1-based, inclusive. We want to make it 0-based, exclusive
    let end: u64 = min(record.end() + padding as u64, contig_len);
    if start > end {
        return Err(BuildError::FetchError(
            record.seqname().to_string(),
            start,
            end,
        ));
    }

    if faidx.fetch(record.seqname(), start, end).is_err() {
        return Err(BuildError::FetchError(
            record.seqname().to_string(),
            start,
            end,
        ));
    } else {
        let mut seq: Vec<u8> = Vec::with_capacity((end - start) as usize);
        let is_rev = match record.strand() {
            Some(strand) => strand.strand_symbol() == "-",
            _ => {
                return Err(BuildError::MissingStrand(
                    record.attributes().get("gene").unwrap().to_string(),
                ))
            }
        };
        match faidx.read(&mut seq) {
            Err(_) => Err(BuildError::FetchError(
                record.seqname().to_string(),
                start,
                end,
            )),
            Ok(_) => {
                if is_rev {
                    Ok(revcomp(seq))
                } else {
                    Ok(seq.to_owned())
                }
            }
        }
    }
}

fn vcf_contig_field(id: &str, length: u64) -> Vec<u8> {
    format!("{}contig=<ID={},length={}>", META, id, length)
        .as_bytes()
        .to_owned()
}

#[cfg(test)]
mod tests {
    use std::io::Write;
    use std::iter::FromIterator;
    use std::str::FromStr;

    use bio::io::fasta::IndexedReader;

    use crate::panel::{Residue, Variant};

    use super::*;

    const PANEL: &str = "tests/cases/panel.tsv";
    const ANNOTATION: &str = "tests/cases/ann.gff3";
    const REF: &str = "tests/cases/ref.fa";

    #[test]
    fn load_annotations_when_no_genes_in_common_returns_empty() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let genes = HashSet::from_iter(vec!["geneX".to_string()]);

        let actual = load_annotations_for_genes(reader, &genes);
        assert!(actual.is_empty())
    }

    #[test]
    fn load_annotations_for_genes_one_gene_in_common() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let genes = HashSet::from_iter(vec!["geneX".to_string(), "dnaA".to_string()]);

        let actual = load_annotations_for_genes(reader, &genes);
        let (actual_gene, actual_record) = actual.iter().next().unwrap();
        assert_eq!(actual_gene, "dnaA");
        assert_eq!(actual_record.end(), &1524u64)
    }

    #[test]
    fn load_annotations_for_genes_is_cds_returns_empty() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tCDS\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let genes = HashSet::from_iter(vec!["geneX".to_string(), "dnaA".to_string()]);

        let actual = load_annotations_for_genes(reader, &genes);
        assert!(actual.is_empty())
    }

    #[test]
    fn extract_gene_from_index_gene_not_in_index() {
        let padding = 0;
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tCDS\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap_err();
        let expected = BuildError::MissingContig("NC_000962.3".to_string());
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_interval_out_of_bounds() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t100\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap_err();
        let expected = BuildError::FetchError("chr1".to_string(), 99, 16);
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_first_base() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t1\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"G";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_too_much_padding_left_wraps_to_start() {
        let padding = 2;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t1\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"GTA";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_too_much_padding_right_wraps_to_end() {
        let padding = 4;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t16\t16\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"ACCCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_padding_start_and_end_exactly_the_same_as_gene() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t16\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"GTAGGCTGAAAACCCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_on_reverse_strand() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t16\t.\t-\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = revcomp(b"GTAGGCTGAAAACCCC");
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_strand() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t16\t.\t.\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap_err();
        let expected = BuildError::MissingStrand("dnaA".to_string());
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_padding_start_same_as_gene_end_minus_one_from_gene_length(
    ) {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t15\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"GTAGGCTGAAAACCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_padding_start_plus_one_from_gene_start_end_same_as_gene(
    ) {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t2\t16\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

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
        let reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let genes = HashSet::from_iter(vec!["dnaA".to_string()]);
        let annotations = load_annotations_for_genes(reader, &genes);
        let tmp = tempfile::NamedTempFile::new().unwrap();

        let header = builder.create_vcf_header(&annotations);

        let writer =
            bcf::Writer::from_path(tmp.path(), &header, false, bcf::Format::VCF)
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
            reference_file: REF.parse().unwrap(),
            padding: 100,
            outdir: PathBuf::from(outdir.path()),
            match_len: 5,
            force: false,
        };
        let result = builder.run();
        assert!(result.is_ok());

        let mut file1 = std::fs::File::open("tests/cases/expected/dr.prg").unwrap();
        let mut file2 =
            std::fs::File::open(format!("{}/dr.prg", outdir.path().to_string_lossy()))
                .unwrap();

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
            reference_file: REF.parse().unwrap(),
            padding: 100,
            outdir: outdir.to_owned(),
            match_len: 5,
            force: false,
        };
        let result = builder.run();
        assert!(result.is_ok());

        let mut file1 = std::fs::File::open("tests/cases/expected/dr.prg").unwrap();
        let mut file2 =
            std::fs::File::open(format!("{}/dr.prg", outdir.to_string_lossy()))
                .unwrap();

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
    fn build_runner_with_force() {
        let outdir = tempfile::tempdir().unwrap();
        let mut builder = Build {
            bcftools_exec: Some(PathBuf::from("src/ext/bcftools")),
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            makeprg_exec: Some(PathBuf::from("src/ext/make_prg")),
            mafft_exec: Some(PathBuf::from("src/ext/mafft/bin/mafft")),
            gff_file: ANNOTATION.parse().unwrap(),
            panel_file: PANEL.parse().unwrap(),
            reference_file: REF.parse().unwrap(),
            padding: 100,
            outdir: PathBuf::from(outdir.path()),
            match_len: 5,
            force: true,
        };
        let result = builder.run();
        assert!(result.is_ok());

        let mut file1 = std::fs::File::open("tests/cases/expected/dr.prg").unwrap();
        let mut file2 =
            std::fs::File::open(format!("{}/dr.prg", outdir.path().to_string_lossy()))
                .unwrap();

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
}
