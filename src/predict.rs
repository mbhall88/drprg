use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::string::ToString;

use anyhow::{Context, Result};
use log::{debug, info};
use rust_htslib::bcf;
use rust_htslib::bcf::header::{TagLength, TagType};
use rust_htslib::bcf::{Format, Read};
use serde::{de, Deserialize, Deserializer, Serialize};
use structopt::clap::AppSettings;
use structopt::StructOpt;
use strum_macros::EnumString;
use thiserror::Error;
use uuid::Uuid;

use drprg::filter::{Filter, Filterer};
use drprg::{unwrap_or_continue, Pandora, PathExt, VcfExt};

use crate::cli::check_path_exists;
use crate::Runner;

/// A collection of custom errors relating to the predict component of this package
#[derive(Error, Debug, PartialEq)]
pub enum PredictError {
    /// The index is not valid
    #[error("Index is not valid due to missing file {0:?}")]
    InvalidIndex(PathBuf),
}

/// All possible predictions
#[derive(Debug, Eq, PartialEq, EnumString, strum_macros::Display, Serialize)]
pub enum Prediction {
    #[strum(to_string = "S")]
    #[serde(alias = "S", rename(serialize = "S"))]
    Susceptible,
    #[strum(to_string = "R")]
    #[serde(alias = "R", rename(serialize = "R"))]
    Resistant,
    #[strum(to_string = "F")]
    #[serde(alias = "F", rename(serialize = "F"))]
    Failed,
    #[strum(to_string = "U")]
    #[serde(alias = "U", rename(serialize = "U"))]
    Unknown,
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
    fn to_bytes(&self) -> &[u8] {
        match self {
            Self::Susceptible => b"S",
            Self::Resistant => b"R",
            Self::Failed => b"F",
            Self::Unknown => b"U",
        }
    }
}

#[derive(StructOpt, Debug, Default)]
#[structopt(setting = AppSettings::DeriveDisplayOrder)]
pub struct Predict {
    /// Path to pandora executable. Will try in src/ext or $PATH if not given
    #[structopt(
        short = "p",
        long = "pandora",
        parse(from_os_str),
        hidden_short_help = true,
        value_name = "FILE"
    )]
    pandora_exec: Option<PathBuf>,
    /// Directory containing the index (produced by `drprg build`)
    #[structopt(short = "x", long, required = true, parse(try_from_os_str = check_path_exists), value_name = "DIR")]
    index: PathBuf,
    /// Sample reads to predict resistance from
    ///
    /// Both fasta and fastq are accepted, along with compressed or uncompressed.
    #[structopt(short, long, required = true, parse(try_from_os_str = check_path_exists), value_name = "FILE")]
    input: PathBuf,
    /// Directory to place output
    #[structopt(
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
    #[structopt(short, long)]
    sample: Option<String>,
    /// Attempt to discover novel variants (i.e. variants not in the panel)
    ///
    /// If a novel variant is discovered, a prediciton of "unknown" is returned for the drug(s)
    /// associated with that gene
    #[structopt(short = "u", long = "unknown")]
    discover: bool,
    /// Require all resistance-associated sites to be genotyped
    ///
    /// If a genotype cannot be assigned for a site (i.e. a null call), then a prediction of
    /// "failed" is returned for the drug(s) associated with that site.
    #[structopt(short = "f", long = "failed")]
    require_genotype: bool,
    /// Sample reads are from Illumina sequencing
    #[structopt(short = "I", long = "illumina")]
    is_illumina: bool,
    #[structopt(flatten)]
    filterer: Filterer,
}

impl Runner for Predict {
    fn run(&mut self) -> Result<()> {
        // todo: implement force
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

        if self.discover {
            todo!("run pandora discover");
        }

        info!("Genotyping reads against the panel with pandora");
        let pandora = Pandora::from_path(&self.pandora_exec)?;
        let threads = &rayon::current_num_threads().to_string();
        let mut gt_args = vec!["-t", threads];
        if self.is_illumina {
            gt_args.push("-I");
        }
        pandora.genotype_with(
            &self.index_prg_path(),
            &self.index_vcf_ref_path(),
            &self.input,
            &self.outdir,
            &gt_args,
        )?;
        info!("Successfully genotyped reads");
        debug!("Filtering variants...");
        let pandora_vcf_path = self.outdir.join(Pandora::vcf_filename());
        let _predict_vcf_path = self.predict_from_pandora_vcf(&pandora_vcf_path)?;
        debug!("Variant filtering complete");

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

    fn index_prg_index_path(&self) -> PathBuf {
        self.index.join("dr.prg.k15.w14.idx")
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

    fn sample_name(&self) -> &str {
        match &self.sample {
            Some(s) => s.as_str(),
            None => self.input.file_prefix().unwrap(),
        }
    }

    fn predict_vcf_filename(&self) -> PathBuf {
        let path = PathBuf::from(self.sample_name());
        path.add_extension(".drprg.bcf".as_ref())
    }

    fn validate_index(&self) -> Result<(), PredictError> {
        let expected_paths = &[
            self.index_prg_path(),
            self.index_kmer_prgs_path(),
            self.index_vcf_path(),
            self.index_vcf_index_path(),
            self.index_vcf_ref_path(),
            self.index_prg_index_path(),
            self.index_prg_update_path(),
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
            bcf::Writer::from_path(&predict_vcf_path, &vcf_header, false, Format::BCF)
                .context("Failed to create filtered VCF")?;
        let mut vcfidx = bcf::IndexedReader::from_path(self.index_vcf_path())
            .context("Failed to open panel VCF index")?;

        for (i, record_result) in reader.records().enumerate() {
            let mut record = record_result
                .context(format!("Failed to read record {} in pandora VCF", i))?;

            writer.translate(&mut record);
            self.filterer.filter(&mut record)?;
            record
                .set_id(&Uuid::new_v4().to_string()[..8].as_bytes())
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
            unwrap_or_continue!(vcfidx.fetch(idx_rid, iv.start as u64, iv.end as u64));
            let mut record_has_resistant = false;
            for idx_res in vcfidx.records() {
                let idx_record = idx_res.context("Failed to parse index vcf record")?;
                let mut prediction = Prediction::Susceptible;
                if record.called_allele() == -1 && self.require_genotype {
                    prediction = Prediction::Failed;
                } else {
                    match record.argmatch(&idx_record) {
                        None if self.discover => {
                            if record_has_resistant {
                                prediction = Prediction::Susceptible
                            } else {
                                prediction = Prediction::Unknown
                            }
                        }
                        Some(i) if i > 0 => {
                            record_has_resistant = true;
                            prediction = Prediction::Resistant
                        }
                        _ => (),
                    };
                }
                // todo: improve this. need to find a way to push all of these into a vector and just have a single update at the end of the loop
                // I've tried, but the borrow checker is too smart for me and I can't solve the problem
                // This works for now, but slows the execution down quite a bit
                let vid = idx_record.id();
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
                        if record_has_resistant {
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
        // todo: generate prediction from filtered pandora output
        Ok(predict_vcf_path)
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

    #[test]
    fn sample_name_no_sample() {
        let predictor = Predict {
            pandora_exec: None,
            index: Default::default(),
            input: PathBuf::from("foo/sample1.fq.gz"),
            outdir: Default::default(),
            sample: None,
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };

        let actual = predictor.sample_name();
        let expected = "sample1";

        assert_eq!(actual, expected)
    }

    #[test]
    fn sample_name_with_sample() {
        let predictor = Predict {
            pandora_exec: None,
            index: Default::default(),
            input: PathBuf::from("foo/sample1.fq.gz"),
            outdir: Default::default(),
            sample: Some("sample2".to_string()),
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
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
            input: Default::default(),
            outdir: Default::default(),
            sample: None,
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
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
            input: Default::default(),
            outdir: Default::default(),
            sample: None,
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };

        let actual = predictor.index_prg_index_path();
        let expected = PathBuf::from("foo/dr.prg.k15.w14.idx");

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_prg_update_path() {
        let predictor = Predict {
            pandora_exec: None,
            sample: None,
            index: PathBuf::from("foo"),
            input: Default::default(),
            outdir: Default::default(),
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };

        let actual = predictor.index_prg_update_path();
        let expected = PathBuf::from("foo/dr.update_DS");

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_kmer_prgs_path() {
        let predictor = Predict {
            pandora_exec: None,
            index: PathBuf::from("foo"),
            input: Default::default(),
            sample: None,
            outdir: Default::default(),
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };

        let actual = predictor.index_kmer_prgs_path();
        let expected = PathBuf::from("foo/kmer_prgs");

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_vcf_path() {
        let predictor = Predict {
            pandora_exec: None,
            index: PathBuf::from("foo"),
            input: Default::default(),
            outdir: Default::default(),
            sample: None,
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };

        let actual = predictor.index_vcf_path();
        let expected = PathBuf::from("foo/panel.bcf");

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_vcf_ref_path() {
        let predictor = Predict {
            pandora_exec: None,
            index: PathBuf::from("foo"),
            input: Default::default(),
            sample: None,
            outdir: Default::default(),
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
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
            pandora_exec: None,
            index: PathBuf::from(tmp_path),
            input: Default::default(),
            outdir: Default::default(),
            sample: None,
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
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
        assert!(predictor.validate_index().is_ok())
    }

    #[test]
    fn validate_index_missing_prg() {
        let predictor = Predict {
            pandora_exec: None,
            index: Default::default(),
            input: Default::default(),
            outdir: Default::default(),
            discover: false,
            sample: None,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(PathBuf::new().join("dr.prg"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_kmer_prgs() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            pandora_exec: None,
            index: PathBuf::from(tmp_path),
            input: Default::default(),
            outdir: Default::default(),
            sample: None,
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
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
            pandora_exec: None,
            index: PathBuf::from(tmp_path),
            sample: None,
            input: Default::default(),
            outdir: Default::default(),
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
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
            pandora_exec: None,
            index: PathBuf::from(tmp_path),
            sample: None,
            input: Default::default(),
            outdir: Default::default(),
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
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
            pandora_exec: None,
            index: PathBuf::from(tmp_path),
            input: Default::default(),
            outdir: Default::default(),
            sample: None,
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf.csi")).unwrap();
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
            pandora_exec: None,
            index: PathBuf::from(tmp_path),
            input: Default::default(),
            outdir: Default::default(),
            discover: false,
            sample: None,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
        };
        {
            let _f = File::create(tmp_path.join("dr.prg")).unwrap();
            let _f = File::create(tmp_path.join("kmer_prgs")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf")).unwrap();
            let _f = File::create(tmp_path.join("panel.bcf.csi")).unwrap();
            let _f = File::create(tmp_path.join("genes.fa")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("dr.prg.k15.w14.idx"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn validate_index_missing_prg_update() {
        let dir = tempfile::tempdir().unwrap();
        let tmp_path = dir.path();
        let predictor = Predict {
            pandora_exec: None,
            index: PathBuf::from(tmp_path),
            input: Default::default(),
            outdir: Default::default(),
            sample: None,
            discover: false,
            require_genotype: false,
            is_illumina: false,
            filterer: Default::default(),
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
    fn predict_add_predict_info_to_header() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        let pred = Predict::default();
        pred.add_predict_info_to_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            ..Default::default()
        };
        let pandora_vcf_path = Path::new("tests/cases/predict/in.bcf");

        let result = pred.predict_from_pandora_vcf(&pandora_vcf_path);
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
            discover: true,
            require_genotype: true,
            ..Default::default()
        };
        let pandora_vcf_path = Path::new("tests/cases/predict/in.bcf");

        let result = pred.predict_from_pandora_vcf(&pandora_vcf_path);
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
}
