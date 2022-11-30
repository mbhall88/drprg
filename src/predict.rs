use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::string::ToString;

use anyhow::{anyhow, Context, Result};
use clap::{AppSettings, Parser};
use log::{debug, info};
use rust_htslib::bcf;
use rust_htslib::bcf::header::{TagLength, TagType};
use rust_htslib::bcf::{Format, Read};
use serde::{de, Deserialize, Deserializer};
use serde_derive::Serialize;
use serde_json::json;
use std::collections::BTreeMap;
use strum_macros::EnumString;
use thiserror::Error;
use uuid::Uuid;

use drprg::filter::{Filter, Filterer};
use drprg::{
    find_prg_index_in, unwrap_or_continue, MakePrg, MultipleSeqAligner, Pandora,
    PathExt, VcfExt,
};

use crate::cli::check_path_exists;
use crate::panel::{Residue, Variant};
use crate::report::{Evidence, Susceptibility};
use crate::Runner;
use bstr::ByteSlice;

use std::collections::{HashMap, HashSet};

use crate::config::Config;
use crate::consequence::consequence_of_variant;
use crate::expert::{ExpertRules, RuleExt};
use crate::minor::MinorAllele;
use noodles::fasta;
use regex::Regex;
use std::fs::File;
use std::io::{BufReader, Write};
use std::iter::FromIterator;

static NONE_DRUG: &str = "NONE";

/// A collection of custom errors relating to the predict component of this package
#[derive(Error, Debug, PartialEq, Eq)]
pub enum PredictError {
    /// The index is not valid
    #[error("Index is not valid due to missing file {0:?}")]
    InvalidIndex(PathBuf),
    /// Could not retrieve the kmer and window size of the index
    #[error("Could not retrieve the kmer and window size of the index PRG")]
    MissingKmerAndWindowSize,
    /// Issue with trying to generate consequence
    #[error("Could not determine consequence of variant because: {0:?}")]
    ConsequenceError(String),
    /// Ran into an error when checking for minor allele
    #[error("{0}")]
    MinorAlleleIssue(String),
}

/// All possible predictions
#[derive(
    Debug,
    Eq,
    PartialEq,
    EnumString,
    strum_macros::Display,
    Serialize,
    Copy,
    Clone,
    PartialOrd,
    Ord,
)]
pub enum Prediction {
    #[strum(to_string = ".")]
    #[serde(alias = ".", rename(serialize = "."))]
    None,
    #[strum(to_string = "S")]
    #[serde(alias = "S", rename(serialize = "S"))]
    Susceptible,
    #[strum(to_string = "F")]
    #[serde(alias = "F", rename(serialize = "F"))]
    Failed,
    #[strum(to_string = "u")]
    #[serde(alias = "u", rename(serialize = "u"))]
    MinorUnknown,
    #[strum(to_string = "U")]
    #[serde(alias = "U", rename(serialize = "U"))]
    Unknown,
    #[strum(to_string = "r")]
    #[serde(alias = "r", rename(serialize = "r"))]
    MinorResistant,
    #[strum(to_string = "R")]
    #[serde(alias = "R", rename(serialize = "R"))]
    Resistant,
}

impl<'de> Deserialize<'de> for Prediction {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        FromStr::from_str(&s).map_err(de::Error::custom)
    }
}

impl Default for Prediction {
    fn default() -> Self {
        Self::None
    }
}

impl Prediction {
    fn to_bytes(self) -> &'static [u8] {
        match self {
            Self::None => b".",
            Self::Susceptible => b"S",
            Self::Resistant => b"R",
            Self::Failed => b"F",
            Self::Unknown => b"U",
            Self::MinorResistant => b"r",
            Self::MinorUnknown => b"u",
        }
    }
}

#[derive(Parser, Debug, Default)]
#[clap(setting = AppSettings::DeriveDisplayOrder)]
pub struct Predict {
    /// Path to pandora executable. Will try in src/ext or $PATH if not given
    #[clap(
        short = 'p',
        long = "pandora",
        parse(from_os_str),
        hidden_short_help = true,
        value_name = "FILE"
    )]
    pandora_exec: Option<PathBuf>,
    /// Path to make_prg executable. Will try in src/ext or $PATH if not given
    #[clap(
        short = 'm',
        long = "makeprg",
        parse(from_os_str),
        hidden_short_help = true,
        value_name = "FILE"
    )]
    makeprg_exec: Option<PathBuf>,
    /// Path to MAFFT executable. Will try in src/ext or $PATH if not given
    #[clap(
        short = 'M',
        long = "mafft",
        parse(from_os_str),
        hidden_short_help = true,
        value_name = "FILE"
    )]
    mafft_exec: Option<PathBuf>,
    /// Directory containing the index (produced by `drprg build`)
    #[clap(short = 'x', long, required = true, parse(try_from_os_str = check_path_exists), value_name = "DIR")]
    index: PathBuf,
    /// Sample reads to predict resistance from
    ///
    /// Both fasta and fastq are accepted, along with compressed or uncompressed.
    #[clap(short, long, required = true, parse(try_from_os_str = check_path_exists), value_name = "FILE")]
    input: PathBuf,
    /// Directory to place output
    #[clap(
        short,
        long,
        default_value = ".",
        parse(from_os_str),
        value_name = "DIR"
    )]
    outdir: PathBuf,
    /// Identifier to use for the sample
    ///
    /// If not provided, this will be set to the input reads file path prefix
    #[clap(short, long)]
    sample: Option<String>,
    /// Sample reads are from Illumina sequencing
    #[clap(short = 'I', long = "illumina")]
    is_illumina: bool,
    /// Minimum allele frequency to call variants
    ///
    /// If an alternate allele has at least this fraction of the depth, a minor resistance
    /// ("r") prediction is made. Set to 1 to disable
    #[clap(
        short = 'f',
        long = "maf",
        value_name = "FLOAT[0.0-1.0]",
        default_value = "1.0"
    )]
    min_allele_freq: f32,
    /// Ignore unknown (off-catalogue) variants that cause a synonymous substitution
    #[clap(short = 'S', long)]
    ignore_synonymous: bool,
    #[clap(flatten)]
    filterer: Filterer,
    /// Set the minimum cluster size in pandora
    #[clap(short = 'C', long, hidden = true, default_value = "10")]
    pandora_min_cluster_size: u16,
    /// Output debugging files. Mostly for development purposes
    #[clap(long, hidden_short_help = true)]
    debug: bool,
}

impl Runner for Predict {
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

        self.validate_index()?;
        let config = Config::from_path(&self.index_config())
            .context("Failed to load index config file")?;
        debug!("Index is valid");
        let threads = &rayon::current_num_threads().to_string();
        let pandora = Pandora::from_path(&self.pandora_exec)?;

        info!("Discovering variants...");
        let tsvpath = self.outdir.join("query.tsv");
        {
            let mut file =
                File::create(&tsvpath).context("Failed to create sample TSV file")?;
            let content = format!(
                "{}\t{}\n",
                self.sample_name(),
                self.input.canonicalize()?.to_string_lossy()
            );
            file.write_all(content.as_bytes())
                .context("Failed to write content to TSV file")?;
        }
        // safe to unwrap as we have already validated the index, which checks this
        let w = config.w.to_string();
        let k = config.k.to_string();
        let c = self.pandora_min_cluster_size.to_string();
        let mut args = vec!["-t", threads, "-w", &w, "-k", &k, "-c", &c];
        if self.is_illumina {
            args.push("-I");
        }
        if self.debug {
            args.push("-K");
        }

        debug!("Running pandora discover...");
        let denovo_paths = pandora
            .discover_with(
                &self.index_prg_path(),
                &tsvpath,
                &self.discover_dir(),
                &args,
            )
            .context("Failed to run pandora discover")?;
        std::fs::remove_file(tsvpath).context("Failed to delete TSV file")?;

        // ==============================================================================================
        debug!("Updating the PRGs with novel variants");
        let makeprg = MakePrg::from_path(&self.makeprg_exec)?;
        let mafft = MultipleSeqAligner::from_path(&self.mafft_exec)?;
        let update_args = &[
            "-t",
            threads,
            "--min_match_len",
            &config.min_match_len.to_string(),
            "--max_nesting",
            &config.max_nesting.to_string(),
        ];
        let prg_path = makeprg
            .update(
                &self.index_msa_dir().canonicalize()?,
                &self.index_prg_path().canonicalize()?,
                &denovo_paths.canonicalize()?,
                &self.outdir.canonicalize()?,
                update_args,
                mafft.executable,
            )
            .context("Failed to update discover PRG")?;
        // ==============================================================================================

        debug!("Indexing updated PRG with pandora index...");
        pandora.index_with(&prg_path, ["-t", threads, "-w", &w, "-k", &k])?;

        let pandora_vcf_path = self.outdir.join(Pandora::vcf_filename());
        info!("Genotyping reads against the panel with pandora");
        // safe to unwrap as we have already validated the index, which checks this
        let mut gt_args = vec!["-t", threads, "-w", &w, "-k", &k, "-c", &c];
        if self.is_illumina {
            gt_args.push("-I");
        }
        if self.debug {
            args.push("-K");
        }

        pandora.genotype_with(
            &prg_path,
            &self.index_vcf_ref_path(),
            &self.input,
            &self.outdir,
            &gt_args,
        )?;
        info!("Successfully genotyped reads");

        info!("Making predictions from variants...");
        let predict_vcf_path = self.predict_from_pandora_vcf(&pandora_vcf_path)?;
        debug!("Predictions written to VCF {:?}", predict_vcf_path);
        let predict_json_path = self.vcf_to_json(&predict_vcf_path)?;
        info!(
            "Prediction report written to JSON file {:?}",
            predict_json_path
        );

        Ok(())
    }
}

impl Predict {
    fn add_predict_info_to_header(&self, header: &mut bcf::Header) {
        for field in &[InfoField::VarId, InfoField::Prediction] {
            let line = format!(
                r#"##INFO=<ID={},Number={},Type={},Description="{}">"#,
                field.id(),
                field.number().as_char(),
                field.tag_type().as_str(),
                field.description()
            );
            header.push_record(line.as_bytes());
        }
    }
    fn index_prg_path(&self) -> PathBuf {
        self.index.join("dr.prg")
    }

    fn index_config(&self) -> PathBuf {
        self.index.join(".config.toml")
    }

    fn index_prg_index_path(&self) -> Result<PathBuf, PredictError> {
        find_prg_index_in(&self.index).ok_or_else(|| {
            PredictError::InvalidIndex(self.index.join("dr.prg.kX.wY.idx"))
        })
    }

    fn load_rules(&self) -> Result<ExpertRules, anyhow::Error> {
        let p = self.index.join("rules.csv");
        if p.exists() {
            ExpertRules::from_csv(&p)
        } else {
            Ok(ExpertRules::new())
        }
    }

    fn index_kmer_prgs_path(&self) -> PathBuf {
        self.index.join("kmer_prgs")
    }

    fn index_vcf_path(&self) -> PathBuf {
        self.index.join("panel.bcf")
    }

    fn index_vcf_index_path(&self) -> PathBuf {
        self.index_vcf_path().add_extension(".csi".as_ref())
    }

    fn index_vcf_ref_path(&self) -> PathBuf {
        self.index.join("genes.fa")
    }

    fn index_msa_dir(&self) -> PathBuf {
        self.index.join("msas")
    }

    fn index_prgs_path(&self) -> PathBuf {
        self.index.join("dr_prgs")
    }

    fn sample_name(&self) -> &str {
        match &self.sample {
            Some(s) => s.as_str(),
            None => self.input.prefix().unwrap(),
        }
    }

    fn predict_vcf_filename(&self) -> PathBuf {
        let path = PathBuf::from(self.sample_name());
        path.add_extension(".drprg.bcf".as_ref())
    }

    fn discover_dir(&self) -> PathBuf {
        self.outdir.join("discover")
    }

    fn json_filename(&self) -> PathBuf {
        let path = PathBuf::from(self.sample_name());
        path.add_extension(".drprg.json".as_ref())
    }

    fn validate_index(&self) -> Result<(), PredictError> {
        let prg_index = self.index_prg_index_path()?;
        let expected_paths = &[
            self.index_config(),
            self.index_prg_path(),
            self.index_kmer_prgs_path(),
            self.index_vcf_path(),
            self.index_vcf_index_path(),
            self.index_vcf_ref_path(),
            prg_index,
            self.index_prgs_path(),
            self.index_msa_dir(),
        ];
        for p in expected_paths {
            if !p.exists() {
                return Err(PredictError::InvalidIndex(p.to_owned()));
            }
        }
        Ok(())
    }

    fn predict_from_pandora_vcf(&self, pandora_vcf_path: &Path) -> Result<PathBuf> {
        let predict_vcf_path = self.outdir.join(self.predict_vcf_filename());
        let mut reader = bcf::Reader::from_path(pandora_vcf_path)
            .context("Failed to open pandora VCF")?;
        let mut vcf_header = bcf::Header::from_template(reader.header());
        self.filterer.add_filter_headers(&mut vcf_header);
        self.add_predict_info_to_header(&mut vcf_header);
        let maf_checker = MinorAllele::new(self.min_allele_freq);
        if self.min_allele_freq < 1.0 {
            maf_checker.add_vcf_headers(&mut vcf_header);
        }
        let mut writer =
            bcf::Writer::from_path(&predict_vcf_path, &vcf_header, false, Format::Bcf)
                .context("Failed to create filtered VCF")?;
        let mut vcfidx = bcf::IndexedReader::from_path(self.index_vcf_path())
            .context("Failed to open panel VCF index")?;
        debug!("Loading the panel...");
        let panel = self.load_var_to_drugs()?;
        let expert_rules = self
            .load_rules()
            .context("Failed to load expert rules in index")?;
        debug!("Loaded panel");

        for (i, record_result) in reader.records().enumerate() {
            let mut record = record_result
                .context(format!("Failed to read record {} in pandora VCF", i))?;

            let mut record_mutations: Vec<String> = Vec::new();
            let mut record_predictions: Vec<Prediction> = Vec::new();

            writer.translate(&mut record);
            self.filterer.filter(&mut record)?;
            record
                .set_id(Uuid::new_v4().to_string()[..8].as_bytes())
                .context("Duplicate ID found - 1/270,000,000 chance of this happening - buy a lottery ticket!")?;
            let chrom = record.contig().into_bytes();
            let has_minor_allele =
                match maf_checker.check_for_minor_alternate(&mut record) {
                    Ok(i) if i > 0 => {
                        MinorAllele::adjust_genotype(&mut record, i as i32).context(
                            format!(
                                "Failed to ajust genotype in record {}:{}",
                                record.contig(),
                                record.pos()
                            ),
                        )?;
                        true
                    }
                    Ok(_) => false,
                    Err(_) => {
                        return Err(anyhow!(
                            "Failed to check for minor allele in record {}:{}",
                            record.contig(),
                            record.pos()
                        ))
                    }
                };
            let idx_rid = unwrap_or_continue!(vcfidx.header().name2rid(&chrom));
            let iv = record.range();
            // fetch start and end are BOTH 0-based INCLUSIVE
            unwrap_or_continue!(vcfidx.fetch(
                idx_rid,
                iv.start as u64,
                Some((iv.end - 1) as u64)
            ));
            let ev = self.consequence(&record).context(format!(
                "Failed to find the consequence of novel variant {}",
                record.id().to_str_lossy()
            ))?;
            let csqs = ev.atomise();

            for idx_res in vcfidx.records() {
                let idx_record = idx_res.context("Failed to parse index vcf record")?;
                let vid = idx_record.id();
                let vid_str = vid.to_str_lossy();
                let (drugs, _) = &panel.get(&*vid_str).unwrap();
                let mut prediction = Prediction::None;
                if record.called_allele() == -1 {
                    prediction = Prediction::Failed;
                } else {
                    match record.argmatch(&idx_record) {
                        None => {
                            for csq in &csqs {
                                let is_x_mutation = vid_str.ends_with('X');
                                let csq_str = csq.to_variant_string();
                                let csq_matches_variant = if is_x_mutation {
                                    csq_str[..csq_str.len() - 1]
                                        == vid_str[..vid_str.len() - 1]
                                } else {
                                    csq_str == vid_str
                                };
                                if csq_matches_variant {
                                    if !drugs.contains(NONE_DRUG) {
                                        prediction = Prediction::Resistant;
                                        break;
                                    } else {
                                        prediction = Prediction::Susceptible;
                                        break;
                                    }
                                }
                            }
                        }
                        Some(i) if i > 0 => {
                            if !drugs.contains(NONE_DRUG) {
                                prediction = Prediction::Resistant
                            }
                        }
                        _ => (),
                    };
                }
                record_predictions.push(prediction);
                record_mutations.push(vid_str.to_string());
            }
            for csq in &csqs {
                let var_str = csq.to_variant_string();
                let mut pred = Prediction::Susceptible;
                let rule_matches = expert_rules.matches(csq);
                if rule_matches.is_empty() {
                    continue;
                }
                for rule in rule_matches {
                    if !rule.drugs.contains(NONE_DRUG) {
                        match record.called_allele() {
                            -1 => pred = Prediction::Failed,
                            i if i > 0 => pred = Prediction::Resistant,
                            _ => pred = Prediction::None,
                        }
                        break;
                    }
                }
                record_mutations.push(var_str);
                record_predictions.push(pred);
            }
            match record_predictions.iter().max() {
                Some(Prediction::None) | None if record.called_allele() > 0 => {
                    for csq in csqs {
                        record_mutations.push(csq.to_variant_string());
                        if csq.is_synonymous() && self.ignore_synonymous {
                            record_predictions.push(Prediction::None);
                        } else {
                            record_predictions.push(Prediction::Unknown);
                        }
                    }
                }
                _ => {}
            }
            if has_minor_allele {
                for p in &mut record_predictions {
                    match p {
                        Prediction::Unknown => *p = Prediction::MinorUnknown,
                        Prediction::Resistant => *p = Prediction::MinorResistant,
                        _ => {}
                    }
                }
            }
            (record_mutations, record_predictions) =
                deduplicate_predictions(record_mutations, record_predictions);
            let mut_bytes: Vec<&[u8]> =
                record_mutations.iter().map(|s| s.as_bytes()).collect();
            record
                .push_info_string(InfoField::VarId.id().as_bytes(), &mut_bytes)
                .context("Failed to push VARIDs to record")?;

            let pred_bytes: Vec<&[u8]> =
                record_predictions.iter().map(|p| p.to_bytes()).collect();
            record
                .push_info_string(InfoField::Prediction.id().as_bytes(), &pred_bytes)
                .context("Failed to push PREDs to record")?;

            writer
                .write(&record)
                .context("Failed to write filtered VCF record")?;
        }
        Ok(predict_vcf_path)
    }

    fn load_var_to_drugs(&self) -> Result<HashMap<String, (HashSet<String>, Residue)>> {
        let mut drug_info = HashMap::new();
        let mut reader = bcf::Reader::from_path(&self.index_vcf_path())
            .context("Failed to open panel VCF")?;
        for record_result in reader.records() {
            let record = unwrap_or_continue!(record_result);
            let drugs: Vec<String> = match record
                .info(b"DRUGS")
                .string()
                .context("Couldn't get DRUGS tag from panel VCF")?
            {
                Some(v) => v
                    .iter()
                    .map(|el| String::from_utf8_lossy(el).to_string())
                    .collect(),
                None => continue,
            };
            let res = match record
                .info(b"RES")
                .string()
                .context("Couldn't get DRUGS tag from panel VCF")?
            {
                Some(v) => Residue::from_str(v.clone()[0].to_str_lossy().as_ref()),
                None => Ok(Residue::Nucleic), // this can't happen
            }?;
            // bcf specs ensure no duplicate IDs
            drug_info.insert(
                String::from_utf8_lossy(&record.id()).into_owned(),
                (HashSet::from_iter(drugs), res),
            );
        }

        Ok(drug_info)
    }

    fn vcf_to_json(&self, vcf_path: &Path) -> Result<PathBuf> {
        let json_path = self.outdir.join(self.json_filename());

        debug!("Loading the panel...");
        let var2drugs = self.load_var_to_drugs()?;
        let mut gene2drugs: HashMap<String, HashSet<String>> = HashMap::new();
        for (var, (drugs, _)) in &var2drugs {
            let (chrom, _) = var
                .split_once('_')
                .context(format!("Couldn't split variant ID {} at underscore", var))?;
            let entry = gene2drugs
                .entry(chrom.to_string())
                .or_insert_with(HashSet::new);
            drugs.iter().for_each(|d| {
                entry.insert(d.to_owned());
            });
        }
        let expert_rules = self
            .load_rules()
            .context("Failed to load expert rules in index")?;
        for (gene, rules) in &expert_rules {
            let entry = gene2drugs
                .entry(gene.to_owned())
                .or_insert_with(HashSet::new);
            rules.iter().for_each(|r| {
                for d in &r.drugs {
                    entry.insert(d.to_owned());
                }
            });
        }
        debug!("Loaded panel");

        let mut json: BTreeMap<String, Susceptibility> = BTreeMap::new();
        let mut reader =
            bcf::Reader::from_path(vcf_path).context("Failed to open predict VCF")?;
        for (i, record_result) in reader.records().enumerate() {
            let record = record_result
                .context(format!("Failed to read record {} in predict VCF", i))?;
            let is_alt = record.called_allele() > 0;
            let preds = match record
                .info(InfoField::Prediction.id().as_bytes())
                .string()
                .context(format!(
                    "Failed to get prediction tag for record number {}",
                    i
                ))? {
                Some(v) => v
                    .iter()
                    .map(|el| Prediction::from_str(el.to_str_lossy().as_ref()).unwrap())
                    .collect::<Vec<Prediction>>(),
                None => vec![],
            };
            if preds.is_empty() && is_alt {
                return Err(anyhow!(
                    "{} tag is unexpectedly empty in VCF",
                    InfoField::Prediction.id()
                ));
            }

            let varids = match record
                .info(InfoField::VarId.id().as_bytes())
                .string()
                .context(format!("Failed to get variant ID for record number {}", i))?
            {
                Some(v) => v
                    .iter()
                    .map(|el| String::from_utf8_lossy(el).to_string())
                    .collect::<Vec<String>>(),
                None => vec![],
            };
            if varids.is_empty() && is_alt {
                return Err(anyhow!(
                    "{} tag is unexpectedly empty in VCF",
                    InfoField::VarId.id()
                ));
            }

            // safe to unwrap as we know preds is not empty
            let max_pred = preds.iter().max().unwrap_or(&Prediction::None);

            // we basically ignore the FILTER column if the variant failed genotyping as we want to
            // output this because it could indicate a deletion or some other event
            let is_failed = *max_pred == Prediction::Failed;
            if (!record.is_pass() && !is_failed) || *max_pred == Prediction::None {
                continue;
            }

            for (prediction, varid) in preds
                .iter()
                .zip(varids.iter())
                .filter(|(p, _)| *p == max_pred)
            {
                let (chrom, var) = varid.split_once('_').context(format!(
                    "Couldn't split variant ID {} at underscore",
                    varid
                ))?;
                let (drugs, residue) = match var2drugs.get(varid) {
                    Some((d, r)) => (d.clone(), r.clone()),
                    None => {
                        // check if in expert rules - find the drug(s) for the rule(s)
                        let ev = self.consequence(&record).context(format!(
                            "Failed to find the consequence of novel variant {}",
                            record.id().to_str_lossy()
                        ))?;
                        let csqs = ev.atomise();
                        let mut residue = None;
                        let mut drugs = HashSet::new();
                        for csq in &csqs {
                            let var_str = csq.to_variant_string();
                            if &var_str == varid {
                                for rule in expert_rules.matches(csq) {
                                    for d in rule.drugs {
                                        let _ = drugs.insert(d);
                                    }
                                }
                                residue = Some(csq.residue.clone());
                                break;
                            }
                        }

                        if drugs.is_empty() {
                            // check if in gene2drugs
                            if let Some(d) = gene2drugs.get(chrom) {
                                for s in d {
                                    drugs.insert(s.to_owned());
                                }
                            };
                        }

                        if residue.is_none() {
                            return Err(anyhow!(
                                "Could not find variant {} in panel or expert rules",
                                varid
                            ));
                        }
                        (drugs, residue.unwrap())
                    }
                };
                let ev = Evidence {
                    variant: Variant::from_str(var)?,
                    gene: chrom.to_owned(),
                    residue: residue.to_owned(),
                    vcfid: String::from_utf8_lossy(&record.id()).into_owned(),
                };
                for drug in drugs.iter().filter(|d| *d != NONE_DRUG) {
                    let entry = json
                        .entry(drug.to_string())
                        .or_insert_with(Susceptibility::default);
                    match (entry.predict, *prediction) {
                        (existing_pred, current_pred)
                            if existing_pred < current_pred =>
                        {
                            entry.predict = *prediction;
                            entry.evidence = vec![ev.to_owned()];
                        }
                        (existing_pred, current_pred)
                            if existing_pred == current_pred =>
                        {
                            entry.evidence.push(ev.to_owned());
                        }
                        _ => {}
                    }
                }
            }
        }
        for (_, (drugs, _)) in var2drugs {
            for d in drugs {
                if d != NONE_DRUG {
                    json.entry(d).or_insert_with(Susceptibility::default);
                }
            }
        }
        let data = json!({"sample": self.sample_name(), "susceptibility": json});
        {
            let mut file =
                File::create(&json_path).context("Failed to create JSON file")?;
            let s = serde_json::to_string_pretty(&data)
                .context("Failed to write to JSON file")?;
            write!(&file, "{}", s)?;
            file.flush()?;
        }
        Ok(json_path)
    }

    fn consequence(&self, record: &bcf::Record) -> Result<Evidence, PredictError> {
        let gene_name = record.contig();

        let mut rdr = File::open(self.index_vcf_ref_path())
            .map(BufReader::new)
            .map(fasta::Reader::new)
            .map_err(|_| {
                PredictError::ConsequenceError(
                    "Could not open VCF reference in index".to_string(),
                )
            })?;

        for result in rdr.records() {
            let gene = result.map_err(|_| {
                PredictError::ConsequenceError(
                    "Could not parse record in gene FASTA".to_string(),
                )
            })?;
            if gene.name() != gene_name {
                continue;
            }
            let description = gene.description().ok_or_else(|| {
                PredictError::ConsequenceError(format!(
                    "Expected gene ({}) record to have padding in description",
                    gene_name
                ))
            })?;
            let re = Regex::new(r"padding=(?P<padding>\d+)").unwrap();
            let caps = re.captures(description).ok_or_else(|| {
                PredictError::ConsequenceError(format!(
                    "Could not get padding from header - {}",
                    gene_name
                ))
            })?;
            let padding_as_str = caps
                .name("padding")
                .ok_or_else(|| {
                    PredictError::ConsequenceError(format!(
                        "Could not get padding from header - {}",
                        gene_name
                    ))
                })?
                .as_str();
            // unwrap here as it would not have been captured if it wasn't a positive int
            let padding = i64::from_str(padding_as_str).unwrap();

            return consequence_of_variant(record, padding, &gene)
                .map_err(PredictError::ConsequenceError);
        }

        Err(PredictError::ConsequenceError(format!(
            "Couldn't find gene {} in index FASTA",
            gene_name
        )))
    }
}

/// Remove duplicate mutations and keep the corresponding highest prediction for the mutation
fn deduplicate_predictions(
    mutations: Vec<String>,
    predictions: Vec<Prediction>,
) -> (Vec<String>, Vec<Prediction>) {
    let mut lookup: HashMap<String, Prediction> = HashMap::new();
    for (var, pred) in mutations.iter().zip(predictions) {
        lookup
            .entry(var.to_string())
            .and_modify(|e| *e = pred.max(*e))
            .or_insert(pred);
    }
    let out_mutations: Vec<String> = lookup.keys().map(|m| m.to_owned()).collect();
    let out_predictions: Vec<Prediction> =
        lookup.values().map(|p| p.to_owned()).collect();
    (out_mutations, out_predictions)
}

enum InfoField {
    VarId,
    Prediction,
}

impl InfoField {
    /// The ID used in the VCF INFO header entry
    fn id(&self) -> &str {
        match self {
            InfoField::Prediction => "PREDICT",
            InfoField::VarId => "VARID",
        }
    }
    /// The INFO Type used for this field
    fn tag_type(&self) -> TagType {
        match self {
            InfoField::Prediction => TagType::String,
            InfoField::VarId => TagType::String,
        }
    }
    /// The INFO Number used for this field
    fn number(&self) -> TagLength {
        match self {
            InfoField::Prediction => TagLength::Variable,
            InfoField::VarId => TagLength::Variable,
        }
    }
    /// The INFO Description used for this field
    fn description(&self) -> &str {
        match self {
            InfoField::Prediction => "The drug resistance prediction(s) for the corresponding VARID(s), where 'R' = resistant, 'S' = susceptible, 'F' = failed, and 'U' = unknown",
            InfoField::VarId => "The identifier for the panel variant(s) the record overlaps with"
        }
    }
}

trait TagLengthExt {
    fn as_char(&self) -> char;
}

impl TagLengthExt for TagLength {
    fn as_char(&self) -> char {
        match self {
            Self::Variable => '.',
            Self::Genotypes => 'G',
            Self::Alleles => 'R',
            Self::AltAlleles => 'A',
            Self::Fixed(i) => std::char::from_u32(*i).unwrap_or('.'),
        }
    }
}

trait TagTypeExt {
    fn as_str(&self) -> &str;
}

impl TagTypeExt for TagType {
    fn as_str(&self) -> &str {
        match self {
            Self::String => "String",
            Self::Flag => "Flag",
            Self::Float => "Float",
            Self::Integer => "Integer",
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;

    use rust_htslib::bcf::{Header, Read};
    use tempfile::{NamedTempFile, TempDir};

    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_prediction_ordering() {
        assert!(Prediction::Susceptible < Prediction::Failed);
        assert!(Prediction::Susceptible < Prediction::Unknown);
        assert!(Prediction::Susceptible < Prediction::Resistant);
        assert!(Prediction::Failed < Prediction::Unknown);
        assert!(Prediction::Failed < Prediction::Resistant);
        assert!(Prediction::Unknown < Prediction::Resistant);

        let v = vec![
            Prediction::Susceptible,
            Prediction::Unknown,
            Prediction::Failed,
        ];
        let m = v.iter().max().unwrap();
        assert_eq!(*m, Prediction::Unknown)
    }

    #[test]
    fn sample_name_no_sample() {
        let predictor = Predict {
            pandora_exec: None,
            input: PathBuf::from("foo/sample1.fq.gz"),
            sample: None,
            is_illumina: false,
            ..Default::default()
        };

        let actual = predictor.sample_name();
        let expected = "sample1";

        assert_eq!(actual, expected)
    }

    #[test]
    fn sample_name_with_sample() {
        let predictor = Predict {
            pandora_exec: None,
            input: PathBuf::from("foo/sample1.fq.gz"),
            sample: Some("sample2".to_string()),
            is_illumina: false,
            ..Default::default()
        };

        let actual = predictor.sample_name();
        let expected = "sample2";

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_prg_path() {
        let predictor = Predict {
            pandora_exec: None,
            index: PathBuf::from("foo"),
            sample: None,
            is_illumina: false,
            ..Default::default()
        };

        let actual = predictor.index_prg_path();
        let expected = PathBuf::from("foo/dr.prg");

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_prg_index_path() {
        let predictor = Predict {
            pandora_exec: None,
            index: PathBuf::from("foo"),
            sample: None,
            is_illumina: false,
            ..Default::default()
        };

        let actual = predictor.index_prg_index_path();
        assert!(actual.is_err())
    }

    #[test]
    fn index_prg_update_path() {
        let predictor = Predict {
            index: PathBuf::from("foo"),
            is_illumina: false,
            ..Default::default()
        };

        let actual = predictor.index_prgs_path();
        let expected = PathBuf::from("foo/dr_prgs");

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_kmer_prgs_path() {
        let predictor = Predict {
            index: PathBuf::from("foo"),
            ..Default::default()
        };

        let actual = predictor.index_kmer_prgs_path();
        let expected = PathBuf::from("foo/kmer_prgs");

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_vcf_path() {
        let predictor = Predict {
            index: PathBuf::from("foo"),
            ..Default::default()
        };

        let actual = predictor.index_vcf_path();
        let expected = PathBuf::from("foo/panel.bcf");

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_vcf_ref_path() {
        let predictor = Predict {
            index: PathBuf::from("foo"),
            ..Default::default()
        };

        let actual = predictor.index_vcf_ref_path();
        let expected = PathBuf::from("foo/genes.fa");

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_is_valid() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join(".config.toml")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf.csi")).unwrap();
            let _f = File::create(tmp_path.join("genes.fa")).unwrap();
            let _f = File::create(tmp_path.join("dr.prg.k15.w14.idx")).unwrap();
            std::fs::create_dir(tmp_path.join("dr_prgs")).unwrap();
            std::fs::create_dir(tmp_path.join("msas")).unwrap();
        }
        assert!(predictor.validate_index().is_ok())
    }

    #[test]
    fn validate_index_missing_config() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };
        {
            let _f = File::create(tmp_path.join("dr.prg.k15.w14.idx")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join(".config.toml"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_prg() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };
        {
            let _f = File::create(tmp_path.join("dr.prg.k15.w14.idx")).unwrap();
            let _f = File::create(tmp_path.join(".config.toml")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("dr.prg"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_kmer_prgs() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("dr.prg.k15.w12.idx")).unwrap();
            let _f = File::create(tmp_path.join(".config.toml")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("kmer_prgs"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_panel_vcf() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };

        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("dr.prg.k15.w12.idx")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join(".config.toml")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("panel.bcf"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_panel_vcf_index() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();

        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };

        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
            let _f = File::create(tmp_path.join(".config.toml")).unwrap();
            let _f = File::create(tmp_path.join("dr.prg.k15.w12.idx")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("panel.bcf.csi"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_vcf_ref() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();

        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };

        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf.csi")).unwrap();
            let _f = File::create(tmp_path.join(".config.toml")).unwrap();
            let _f = File::create(tmp_path.join("dr.prg.k15.w12.idx")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("genes.fa"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_prg_index() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf.csi")).unwrap();
            let _f = File::create(tmp_path.join(".config.toml")).unwrap();
            let _f = File::create(tmp_path.join("genes.fa")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("dr.prg.kX.wY.idx"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_msas() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf.csi")).unwrap();
            let _f = File::create(tmp_path.join("genes.fa")).unwrap();
            let _f = File::create(tmp_path.join(".config.toml")).unwrap();
            let _f = File::create(tmp_path.join("dr.prg.k15.w14.idx")).unwrap();
            std::fs::create_dir(tmp_path.join("dr_prgs")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("msas"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_prgs() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            index: PathBuf::from(tmp_path),
            ..Default::default()
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf.csi")).unwrap();
            let _f = File::create(tmp_path.join("genes.fa")).unwrap();
            let _f = File::create(tmp_path.join(".config.toml")).unwrap();
            let _f = File::create(tmp_path.join("dr.prg.k15.w14.idx")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("dr_prgs"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn predict_add_predict_info_to_header() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        let pred = Predict::default();
        pred.add_predict_info_to_header(&mut header);
        let vcf = bcf::Writer::from_path(path, &header, true, Format::Vcf).unwrap();
        let view = vcf.header();

        for field in &[InfoField::VarId, InfoField::Prediction] {
            let result = view.info_type(field.id().as_bytes());
            assert!(result.is_ok());
            let (tag_type, tag_len) = result.unwrap();
            assert_eq!(tag_type, field.tag_type());
            assert_eq!(tag_len, field.number())
        }
    }

    #[test]
    #[allow(clippy::assertions_on_constants)]
    fn test_predict_from_pandora_vcf() {
        let tmp = TempDir::new().unwrap();
        let tmpoutdir = tmp.path();
        let filt = Filterer {
            min_frs: 0.51,
            min_covg: 3,
            min_strand_bias: 0.01,
            max_indel: Some(20),
            min_gt_conf: 5.0,
            ..Default::default()
        };
        let pred = Predict {
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            index: PathBuf::from("tests/cases/predict"),
            outdir: PathBuf::from(tmpoutdir),
            sample: Some("test".to_string()),
            ignore_synonymous: true,
            filterer: filt,
            min_allele_freq: 0.25,
            ..Default::default()
        };
        let pandora_vcf_path = Path::new("tests/cases/predict/in.vcf");

        let result = pred.predict_from_pandora_vcf(pandora_vcf_path);
        assert!(result.is_ok());

        let mut expected_rdr =
            bcf::Reader::from_path(Path::new("tests/cases/predict/out.vcf")).unwrap();
        let mut actual_rdr =
            bcf::Reader::from_path(tmpoutdir.join("test.drprg.bcf")).unwrap();
        let mut actual_records = actual_rdr.records();
        for r in expected_rdr.records() {
            let expected_record = r.unwrap();
            let actual_record = actual_records.next().unwrap().unwrap();
            assert_eq!(expected_record.pos(), actual_record.pos());
            let expected_varid = expected_record.info(b"VARID").string().unwrap();
            let actual_varid = actual_record.info(b"VARID").string().unwrap();
            match (expected_varid, actual_varid) {
                (Some(bb1), Some(bb2)) => {
                    let mut left = bb1.to_owned();
                    let mut right = bb2.to_owned();
                    left.sort_unstable();
                    right.sort_unstable();
                    assert_eq!(
                        left,
                        right,
                        "{}:{} {}:{}",
                        actual_record.contig(),
                        actual_record.pos(),
                        expected_record.contig(),
                        expected_record.pos()
                    )
                }
                (None, Some(_)) | (Some(_), None) => assert!(
                    false,
                    "{}:{} {}:{}",
                    actual_record.contig(),
                    actual_record.pos(),
                    expected_record.contig(),
                    expected_record.pos()
                ),
                (None, None) => assert!(true),
            }

            let expected_pred = expected_record.info(b"PREDICT").string().unwrap();
            let actual_pred = actual_record.info(b"PREDICT").string().unwrap();
            match (expected_pred, actual_pred) {
                (Some(bb1), Some(bb2)) => {
                    let mut left = bb1.to_owned();
                    let mut right = bb2.to_owned();
                    left.sort_unstable();
                    right.sort_unstable();
                    assert_eq!(
                        left,
                        right,
                        "{}:{} {}:{}",
                        actual_record.contig(),
                        actual_record.pos(),
                        expected_record.contig(),
                        expected_record.pos()
                    )
                }
                (None, Some(_)) | (Some(_), None) => assert!(false),
                (None, None) => assert!(true),
            }
        }
    }

    #[test]
    fn test_vcf_to_json() {
        use std::io::Read;
        use std::iter::Iterator;

        let tmp = TempDir::new().unwrap();
        let tmpoutdir = tmp.path();
        let filt = Filterer {
            min_frs: 0.51,
            min_covg: 3,
            min_strand_bias: 0.01,
            max_indel: Some(20),
            min_gt_conf: 5.0,
            ..Default::default()
        };
        let pred = Predict {
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            index: PathBuf::from("tests/cases/predict"),
            outdir: PathBuf::from(tmpoutdir),
            sample: Some("test".to_string()),
            ignore_synonymous: true,
            filterer: filt,
            ..Default::default()
        };
        let vcf_path = Path::new("tests/cases/predict/out.vcf");
        let result = pred.vcf_to_json(vcf_path);
        assert!(result.is_ok());

        let json_path = result.unwrap();
        let file = File::open(json_path).unwrap();
        let mut reader = BufReader::new(file);
        let mut buffer = String::new();
        reader.read_to_string(&mut buffer).unwrap();
        let actual = buffer
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<String>();

        buffer.clear();
        let file = File::open("tests/cases/predict/expected.json").unwrap();
        let mut reader = BufReader::new(file);
        reader.read_to_string(&mut buffer).unwrap();
        let expected = buffer
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<String>();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_deduplicate_predictions() {
        let mutations = vec![
            "embA_V206M".to_string(),
            "embA_V206M".to_string(),
            "embA_V206N".to_string(),
        ];
        let predictions = vec![
            Prediction::Susceptible,
            Prediction::MinorResistant,
            Prediction::Unknown,
        ];

        let (mut actual_mutations, mut actual_predictions) =
            deduplicate_predictions(mutations, predictions);
        let mut expected_mutations =
            vec!["embA_V206M".to_string(), "embA_V206N".to_string()];
        let mut expected_predictions =
            vec![Prediction::MinorResistant, Prediction::Unknown];

        actual_predictions.sort_unstable();
        actual_mutations.sort_unstable();
        expected_mutations.sort_unstable();
        expected_predictions.sort_unstable();

        assert_eq!(actual_predictions, expected_predictions);
        assert_eq!(actual_mutations, expected_mutations)
    }

    #[test]
    fn test_deduplicate_predictions_empty_in_empty_out() {
        let mutations = vec![];
        let predictions = vec![];

        let (actual_mutations, actual_predictions) =
            deduplicate_predictions(mutations, predictions);

        assert!(actual_mutations.is_empty());
        assert!(actual_predictions.is_empty())
    }
}
