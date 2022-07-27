use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::string::ToString;

use anyhow::{Context, Result};
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
    extract_w_and_k, find_prg_index_in, unwrap_or_continue, MakePrg,
    MultipleSeqAligner, Pandora, PathExt, VcfExt,
};

use crate::cli::check_path_exists;
use crate::panel::{Residue, Variant};
use crate::report::{Evidence, Susceptibility};
use crate::Runner;
use bstr::ByteSlice;

use std::collections::{HashMap, HashSet};

use crate::consequence::consequence_of_variant;
use noodles::fasta;
use regex::Regex;
use std::fs::File;
use std::io::{BufReader, Write};
use std::iter::FromIterator;

static NONE_DRUG: &str = "NONE";

/// A collection of custom errors relating to the predict component of this package
#[derive(Error, Debug, PartialEq)]
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
    #[strum(to_string = "S")]
    #[serde(alias = "S", rename(serialize = "S"))]
    Susceptible,
    #[strum(to_string = "F")]
    #[serde(alias = "F", rename(serialize = "F"))]
    Failed,
    #[strum(to_string = "U")]
    #[serde(alias = "U", rename(serialize = "U"))]
    Unknown,
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
        Self::Susceptible
    }
}

impl Prediction {
    fn to_bytes(self) -> &'static [u8] {
        match self {
            Self::Susceptible => b"S",
            Self::Resistant => b"R",
            Self::Failed => b"F",
            Self::Unknown => b"U",
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
    /// Do not output unknown (off-catalogue) variants
    ///
    /// If a novel variant is discovered, a prediction of "unknown" ('U') is returned for the drug(s)
    /// associated with that gene. Using this flag will turn this off an ignore those variants.
    #[clap(short = 'U', long)]
    no_unknown: bool,
    /// Require all resistance-associated sites to be genotyped
    ///
    /// If a genotype cannot be assigned for a site (i.e. a null call), then a prediction of
    /// "failed" ('F') is returned for the drug(s) associated with that site.
    #[clap(short = 'f', long = "failed")]
    require_genotype: bool,
    /// Sample reads are from Illumina sequencing
    #[clap(short = 'I', long = "illumina")]
    is_illumina: bool,
    /// Ignore unknown (off-catalogue) variants that cause a synonymous substitution
    #[clap(short = 'S', long)]
    ignore_synonymous: bool,
    #[clap(flatten)]
    filterer: Filterer,
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
        let (w, k) = self.index_w_and_k().unwrap();
        let mut args = vec!["-t", threads, "-w", &w, "-k", &k];
        if self.is_illumina {
            args.push("-I");
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

        debug!("Running make_prg update...");
        let makeprg = MakePrg::from_path(&self.makeprg_exec)?;
        let mafft = MultipleSeqAligner::from_path(&self.mafft_exec)?;
        let prg_path = makeprg
            .update(
                &self.index_prg_update_path().canonicalize()?,
                &denovo_paths.canonicalize()?,
                &self.outdir,
                &["-t", threads],
                mafft.executable,
            )
            .context("Failed to update discover PRG")?;

        debug!("Indexing updated PRG with pandora index...");
        pandora.index_with(&prg_path, &["-t", threads, "-w", &w, "-k", &k])?;

        let pandora_vcf_path = self.outdir.join(Pandora::vcf_filename());
        info!("Genotyping reads against the panel with pandora");
        // safe to unwrap as we have already validated the index, which checks this
        let (w, k) = self.index_w_and_k().unwrap();
        let mut gt_args = vec!["-t", threads, "-w", &w, "-k", &k];
        if self.is_illumina {
            gt_args.push("-I");
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

    fn index_w_and_k(&self) -> Result<(String, String), PredictError> {
        let fname = self
            .index_prg_index_path()
            .map_err(|_| PredictError::MissingKmerAndWindowSize)?
            .to_string_lossy()
            .to_string();
        let (w, k) =
            extract_w_and_k(&fname).ok_or(PredictError::MissingKmerAndWindowSize)?;
        Ok((w, k))
    }

    fn index_prg_index_path(&self) -> Result<PathBuf, PredictError> {
        find_prg_index_in(&self.index).ok_or_else(|| {
            PredictError::InvalidIndex(self.index.join("dr.prg.kX.wY.idx"))
        })
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

    fn index_prg_update_path(&self) -> PathBuf {
        self.index.join("dr.update_DS")
    }

    fn index_update_prgs_path(&self) -> PathBuf {
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
            self.index_prg_path(),
            self.index_kmer_prgs_path(),
            self.index_vcf_path(),
            self.index_vcf_index_path(),
            self.index_vcf_ref_path(),
            prg_index,
            self.index_prg_update_path(),
            self.index_update_prgs_path(),
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
        let mut writer =
            bcf::Writer::from_path(&predict_vcf_path, &vcf_header, false, Format::Bcf)
                .context("Failed to create filtered VCF")?;
        let mut vcfidx = bcf::IndexedReader::from_path(self.index_vcf_path())
            .context("Failed to open panel VCF index")?;
        debug!("Loading the panel...");
        let panel = self.load_var_to_drugs()?;
        debug!("Loaded panel");

        for (i, record_result) in reader.records().enumerate() {
            let mut record = record_result
                .context(format!("Failed to read record {} in pandora VCF", i))?;

            writer.translate(&mut record);
            self.filterer.filter(&mut record)?;
            record
                .set_id(Uuid::new_v4().to_string()[..8].as_bytes())
                .context("Duplicate ID found - 1/270,000,000 chance of this happening - buy a lottery ticket!")?;
            let rid = record.rid().context(format!(
                "Pandora variant number {} does not have a CHROM",
                i
            ))?;
            let chrom = record
                .header()
                .rid2name(rid)
                .context("Pandora VCF is missing a contig from the header")?;
            let idx_rid = unwrap_or_continue!(vcfidx.header().name2rid(chrom));
            let iv = record.range();
            // fetch start and end are BOTH 0-based INCLUSIVE
            unwrap_or_continue!(vcfidx.fetch(
                idx_rid,
                iv.start as u64,
                Some((iv.end - 1) as u64)
            ));
            let mut record_has_match_in_idx = false;
            for idx_res in vcfidx.records() {
                let idx_record = idx_res.context("Failed to parse index vcf record")?;
                let vid = idx_record.id();
                let mut prediction = Prediction::Susceptible;
                if record.called_allele() == -1 && self.require_genotype {
                    prediction = Prediction::Failed;
                } else {
                    match record.argmatch(&idx_record) {
                        None if !self.no_unknown => {
                            if record_has_match_in_idx {
                                prediction = Prediction::Susceptible
                            } else {
                                prediction = Prediction::Unknown
                            }
                        }
                        Some(i) if i > 0 => {
                            record_has_match_in_idx = true;
                            let vid_str = vid.to_str_lossy();
                            // safe to unwrap here as the var id must be in the panel by definition
                            let (drugs, _) = &panel.get(&*vid_str).unwrap();
                            if !drugs.contains(NONE_DRUG) {
                                prediction = Prediction::Resistant
                            }
                        }
                        _ => (),
                    };
                }
                // todo: improve this. need to find a way to push all of these into a vector and just have a single update at the end of the loop
                // I've tried, but the borrow checker is too smart for me and I can't solve the problem
                // This works for now, but slows the execution down quite a bit
                let pred = prediction.to_bytes();
                let current_varids = record
                    .info(InfoField::VarId.id().as_bytes())
                    .string()
                    .context("Couldn't unwrap VARIDs")?;
                match current_varids {
                    Some(v) => {
                        let mut cur_v = v.clone();
                        cur_v.push(&vid);
                        record
                            .push_info_string(InfoField::VarId.id().as_bytes(), &cur_v)
                            .context("Failed to push VARID to record")?;
                    }
                    None => record
                        .push_info_string(InfoField::VarId.id().as_bytes(), &[&vid])
                        .context("Failed to push VARID to record")?,
                };
                let current_preds = record
                    .info(InfoField::Prediction.id().as_bytes())
                    .string()
                    .context("Couldn't unwrap PREDICTs")?;
                match current_preds {
                    Some(v) => {
                        let mut cur_v = v.clone();
                        if record_has_match_in_idx {
                            // change all unknowns to susceptible
                            let u = Prediction::Unknown.to_bytes();
                            for p in &mut cur_v {
                                if p == &u {
                                    *p = Prediction::Susceptible.to_bytes();
                                }
                            }
                        }
                        cur_v.push(pred);
                        record
                            .push_info_string(
                                InfoField::Prediction.id().as_bytes(),
                                &cur_v,
                            )
                            .context("Failed to push PREDICT to record")?;
                    }
                    None => record
                        .push_info_string(
                            InfoField::Prediction.id().as_bytes(),
                            &[pred],
                        )
                        .context("Failed to push PREDICT to record")?,
                }
            }

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
        let mut var2drugs = self.load_var_to_drugs()?;
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
        debug!("Loaded panel");

        let mut json: BTreeMap<String, Susceptibility> = BTreeMap::new();
        let mut reader =
            bcf::Reader::from_path(vcf_path).context("Failed to open predict VCF")?;
        for (i, record_result) in reader.records().enumerate() {
            let record = record_result
                .context(format!("Failed to read record {} in predict VCF", i))?;
            if !record.is_pass() {
                continue;
            }
            let mut preds = match record
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
            let mut varids = match record
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
            // todo: if the "highest" prediction is unknown, i.e., the variant overlaps stuff in the panel, but doesn't match any, then we clear this list and turn the variant into evidence in the same way we would if the variant didn't overlap the catalogue at all
            let max_pred = preds.iter().max();
            let has_unknown =
                preds.is_empty() || max_pred == Some(&Prediction::Unknown);

            if has_unknown && !self.no_unknown && record.called_allele() > 0 {
                let ev = self.consequence(&record).context(format!(
                    "Failed to find consequence of novel variant {}",
                    record.id().to_str_lossy()
                ))?;

                if !preds.is_empty() {
                    preds.iter_mut().for_each(|p| {
                        if *p == Prediction::Unknown {
                            *p = Prediction::Susceptible
                        } else {
                        }
                    });
                }
                preds.push(Prediction::Unknown);
                let var = format!("{}_{}", ev.gene, ev.variant);
                varids.push(var.to_owned());

                if self.ignore_synonymous && ev.is_synonymous() {
                    continue;
                }

                // for novel variants, add 'U' for all drugs this gene is associated with
                let chrom = record.contig();
                match gene2drugs.get(&chrom) {
                    Some(d) => {
                        if !var2drugs.contains_key(&var) {
                            var2drugs.insert(var, (d.clone(), ev.residue.clone()));
                        }
                        for drug in d.iter().filter(|el| *el != NONE_DRUG) {
                            let entry = json
                                .entry(drug.to_string())
                                .or_insert_with(Susceptibility::default);
                            match entry.predict {
                                Prediction::Susceptible | Prediction::Failed => {
                                    entry.evidence = vec![ev.to_owned()];
                                    entry.predict = Prediction::Unknown;
                                }
                                Prediction::Unknown => {
                                    entry.evidence.push(ev.to_owned())
                                }
                                Prediction::Resistant => {}
                            }
                        }
                    }
                    None => continue,
                }
                continue;
            }
            for (prediction, v) in preds.iter().zip(varids.iter()) {
                if *prediction != Prediction::Susceptible {
                    let (drugs, residue) = var2drugs
                        .get(v)
                        .context(format!("Variant {} in VCF is not in the panel", v))?;
                    let (chrom, var) = v.split_once('_').context(format!(
                        "Couldn't split variant ID {} at underscore",
                        v
                    ))?;
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
                            (p1, p2) if p1 == p2 => entry.evidence.push(ev.to_owned()),
                            (Prediction::Resistant, _) => {}
                            (_, Prediction::Resistant) => {
                                entry.predict = Prediction::Resistant;
                                entry.evidence = vec![ev.to_owned()];
                            }
                            (Prediction::Unknown, _) => {}
                            (_, Prediction::Unknown) => {
                                entry.predict = Prediction::Unknown;
                                entry.evidence = vec![ev.to_owned()];
                            }
                            (Prediction::Failed, _) => {}
                            (_, Prediction::Failed) => {
                                entry.predict = Prediction::Failed;
                                entry.evidence = vec![ev.to_owned()];
                            }
                            _ => {}
                        }
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
            no_unknown: false,
            require_genotype: false,
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
            no_unknown: false,
            require_genotype: false,
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
            no_unknown: false,
            require_genotype: false,
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
            no_unknown: false,
            require_genotype: false,
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
            no_unknown: false,
            require_genotype: false,
            is_illumina: false,
            ..Default::default()
        };

        let actual = predictor.index_prg_update_path();
        let expected = PathBuf::from("foo/dr.update_DS");

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
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf.csi")).unwrap();
            let _f = File::create(tmp_path.join("genes.fa")).unwrap();
            let _f = File::create(tmp_path.join("dr.prg.k15.w14.idx")).unwrap();
            let _f = File::create(tmp_path.join("dr.update_DS")).unwrap();
            std::fs::create_dir(tmp_path.join("dr_prgs")).unwrap();
        }
        assert!(predictor.validate_index().is_ok())
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
            let _f = File::create(tmp_path.join("genes.fa")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("dr.prg.kX.wY.idx"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_prg_update() {
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
            let _f = File::create(tmp_path.join("dr.prg.k15.w14.idx")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("dr.update_DS"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_update_prgs() {
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
            let _f = File::create(tmp_path.join("dr.prg.k15.w14.idx")).unwrap();
            let _f = File::create(tmp_path.join("dr.update_DS")).unwrap();
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
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
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
    fn predict_from_pandora_vcf_no_unknown_or_failed() {
        let tmp = TempDir::new().unwrap();
        let tmpoutdir = tmp.path();
        let filt = Filterer {
            min_frs: 0.7,
            ..Default::default()
        };
        let pred = Predict {
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            index: PathBuf::from("tests/cases/predict"),
            outdir: PathBuf::from(tmpoutdir),
            sample: Some("test".to_string()),
            filterer: filt,
            no_unknown: true,
            ..Default::default()
        };
        let pandora_vcf_path = Path::new("tests/cases/predict/in.bcf");

        let result = pred.predict_from_pandora_vcf(pandora_vcf_path);
        assert!(result.is_ok());

        let mut expected_rdr =
            bcf::Reader::from_path(Path::new("tests/cases/predict/out.bcf")).unwrap();
        let mut actual_rdr =
            bcf::Reader::from_path(tmpoutdir.join("test.drprg.bcf")).unwrap();
        let mut actual_records = actual_rdr.records();
        for r in expected_rdr.records() {
            let expected_record = r.unwrap();
            let actual_record = actual_records.next().unwrap().unwrap();
            assert_eq!(expected_record.pos(), actual_record.pos());
            assert_eq!(
                expected_record
                    .info(b"VARID")
                    .string()
                    .unwrap()
                    .unwrap()
                    .as_slice(),
                actual_record
                    .info(b"VARID")
                    .string()
                    .unwrap()
                    .unwrap()
                    .as_slice()
            );
            assert_eq!(
                expected_record
                    .info(b"PREDICT")
                    .string()
                    .unwrap()
                    .unwrap()
                    .as_slice(),
                actual_record
                    .info(b"PREDICT")
                    .string()
                    .unwrap()
                    .unwrap()
                    .as_slice()
            );
            assert!(expected_record
                .filters()
                .zip(actual_record.filters())
                .all(|(a, b)| expected_record.header().id_to_name(a)
                    == actual_record.header().id_to_name(b)));
        }
    }

    #[test]
    fn predict_from_pandora_vcf_with_unknown_and_failed() {
        let tmp = TempDir::new().unwrap();
        let tmpoutdir = tmp.path();
        let filt = Filterer {
            min_frs: 0.7,
            ..Default::default()
        };
        let pred = Predict {
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            index: PathBuf::from("tests/cases/predict"),
            outdir: PathBuf::from(tmpoutdir),
            sample: Some("test".to_string()),
            filterer: filt,
            no_unknown: false,
            require_genotype: true,
            ..Default::default()
        };
        let pandora_vcf_path = Path::new("tests/cases/predict/in.bcf");

        let result = pred.predict_from_pandora_vcf(pandora_vcf_path);
        assert!(result.is_ok());

        let mut expected_rdr =
            bcf::Reader::from_path(Path::new("tests/cases/predict/out2.bcf")).unwrap();
        let mut actual_rdr =
            bcf::Reader::from_path(tmpoutdir.join("test.drprg.bcf")).unwrap();
        let mut actual_records = actual_rdr.records();
        for r in expected_rdr.records() {
            let expected_record = r.unwrap();
            let actual_record = actual_records.next().unwrap().unwrap();
            assert_eq!(expected_record.pos(), actual_record.pos());
            assert_eq!(
                expected_record
                    .info(b"VARID")
                    .string()
                    .unwrap()
                    .unwrap()
                    .as_slice(),
                actual_record
                    .info(b"VARID")
                    .string()
                    .unwrap()
                    .unwrap()
                    .as_slice()
            );
            assert_eq!(
                expected_record
                    .info(b"PREDICT")
                    .string()
                    .unwrap()
                    .unwrap()
                    .as_slice(),
                actual_record
                    .info(b"PREDICT")
                    .string()
                    .unwrap()
                    .unwrap()
                    .as_slice()
            );
            assert!(expected_record
                .filters()
                .zip(actual_record.filters())
                .all(|(a, b)| expected_record.header().id_to_name(a)
                    == actual_record.header().id_to_name(b)));
        }
    }

    #[test]
    fn vcf_to_json_no_unknown_or_failed() {
        use std::io::Read;
        use std::iter::Iterator;

        let tmp = TempDir::new().unwrap();
        let tmpoutdir = tmp.path();
        let filt = Filterer {
            min_frs: 0.7,
            ..Default::default()
        };
        let pred = Predict {
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            index: PathBuf::from("tests/cases/predict"),
            outdir: PathBuf::from(tmpoutdir),
            sample: Some("test".to_string()),
            filterer: filt,
            ..Default::default()
        };
        let vcf_path = Path::new("tests/cases/predict/out.bcf");
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
        let expected = r#"
        {
          "sample": "test",
          "susceptibility": {
            "Isoniazid": {
              "evidence": [
                {
                  "gene": "fabG1",
                  "residue": "DNA",
                  "variant": "CTG607TTA",
                  "vcfid": "."
                },
                {
                  "gene": "katG",
                  "residue": "DNA",
                  "variant": "G2097GT",
                  "vcfid": "."
                }
              ],
              "predict": "R"
            },
            "Ofloxacin": {
              "evidence": [],
              "predict": "S"
            },
            "Rifampicin": {
              "evidence": [],
              "predict": "S"
            },
            "Streptomycin": {
              "evidence": [
                {
                  "gene": "gid",
                  "residue": "PROT",
                  "variant": "R47W",
                  "vcfid": "."
                },
                {
                  "gene": "rpsL",
                  "residue": "PROT",
                  "variant": "K43R",
                  "vcfid": "."
                }
               ],
              "predict": "R"
            }
          }
        }
        "#
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect::<String>();

        assert_eq!(actual, expected)
    }

    #[test]
    fn vcf_to_json_with_unknown_and_failed() {
        use std::io::Read;
        use std::iter::Iterator;

        let tmp = TempDir::new().unwrap();
        let tmpoutdir = tmp.path();
        let filt = Filterer {
            min_frs: 0.7,
            ..Default::default()
        };
        let pred = Predict {
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            index: PathBuf::from("tests/cases/predict"),
            outdir: PathBuf::from(tmpoutdir),
            sample: Some("test".to_string()),
            filterer: filt,
            no_unknown: true,
            require_genotype: true,
            ..Default::default()
        };
        let vcf_path = Path::new("tests/cases/predict/out2.bcf");
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
        let expected = r#"
        {
          "sample": "test",
          "susceptibility": {
            "Isoniazid": {
              "evidence": [
                {
                  "gene": "fabG1",
                  "residue": "DNA",
                  "variant": "CTG607TTA",
                  "vcfid": "."
                },
                {
                  "gene": "katG",
                  "residue": "DNA",
                  "variant": "G2097GT",
                  "vcfid": "."
                }
              ],
              "predict": "R"
            },
            "Ofloxacin": {
              "evidence": [
                {
                  "gene": "inhA",
                  "residue": "PROT",
                  "variant": "I194T",
                  "vcfid": "."
                }
              ],
              "predict": "F"
            },
            "Rifampicin": {
              "evidence": [
                {
                  "gene": "inhA",
                  "residue": "PROT",
                  "variant": "I21T",
                  "vcfid": "."
                }
              ],
              "predict": "U"
            },
            "Streptomycin": {
              "evidence": [
                {
                  "gene": "gid",
                  "residue": "PROT",
                  "variant": "R47W",
                  "vcfid": "."
                },
                {
                  "gene": "rpsL",
                  "residue": "PROT",
                  "variant": "K43R",
                  "vcfid": "."
                }
               ],
              "predict": "R"
            }
          }
        }
        "#
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect::<String>();

        assert_eq!(actual, expected)
    }

    #[test]
    /// Also, there is a mutation that doesn't overlap the panel that gets added as unknown
    fn vcf_to_json_strep_has_resistant_and_unknown_only_evidence_for_resistant() {
        use std::io::Read;
        use std::iter::Iterator;

        let tmp = TempDir::new().unwrap();
        let tmpoutdir = tmp.path();
        let filt = Filterer {
            min_frs: 0.7,
            ..Default::default()
        };
        let pred = Predict {
            pandora_exec: Some(PathBuf::from("src/ext/pandora")),
            index: PathBuf::from("tests/cases/predict"),
            outdir: PathBuf::from(tmpoutdir),
            sample: Some("test".to_string()),
            filterer: filt,
            no_unknown: false,
            require_genotype: true,
            ..Default::default()
        };
        let vcf_path = Path::new("tests/cases/predict/out3.bcf");
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
        let expected = r#"
        {
          "sample": "test",
          "susceptibility": {
            "Isoniazid": {
              "evidence": [
                {
                  "gene": "fabG1",
                  "residue": "DNA",
                  "variant": "CTG607TTA",
                  "vcfid": "."
                },
                {
                  "gene": "katG",
                  "residue": "DNA",
                  "variant": "G2097GT",
                  "vcfid": "."
                }
              ],
              "predict": "R"
            },
            "Ofloxacin": {
              "evidence": [
                {
                  "gene": "inhA",
                  "residue": "DNA",
                  "variant": "TC62G",
                  "vcfid": "."
                },
                {
                  "gene": "inhA",
                  "residue": "PROT",
                  "variant": "PI227HF",
                  "vcfid": "."
                }
              ],
              "predict": "U"
            },
            "Rifampicin": {
              "evidence": [
                {
                  "gene": "inhA",
                  "residue": "DNA",
                  "variant": "TC62G",
                  "vcfid": "."
                },
                {
                  "gene": "inhA",
                  "residue": "PROT",
                  "variant": "PI227HF",
                  "vcfid": "."
                }
              ],
              "predict": "U"
            },
            "Streptomycin": {
              "evidence": [
                {
                  "gene": "gid",
                  "residue": "PROT",
                  "variant": "R47W",
                  "vcfid": "."
                },
                {
                  "gene": "rpsL",
                  "residue": "PROT",
                  "variant": "K43R",
                  "vcfid": "."
                }
               ],
              "predict": "R"
            }
          }
        }
        "#
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect::<String>();

        assert_eq!(actual, expected)
    }
}
