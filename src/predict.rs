use std::path::PathBuf;

use anyhow::{Context, Result};
use log::{debug, info};
use structopt::clap::AppSettings;
use structopt::StructOpt;
use thiserror::Error;

use crate::cli::check_path_exists;
use crate::Runner;
use drprg::filter::Filterer;
use drprg::{Pandora, PathExt};
use rust_htslib::bcf::Read;

/// A collection of custom errors relating to the predict component of this package
#[derive(Error, Debug, PartialEq)]
pub enum PredictError {
    /// The index is not valid
    #[error("Index is not valid due to missing file {0:?}")]
    InvalidIndex(PathBuf),
}

#[derive(StructOpt, Debug)]
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
    fn run(&self) -> Result<()> {
        if !self.outdir.exists() {
            info!("Outdir doesn't exist...creating...");
            std::fs::create_dir(&self.outdir)
                .context(format!("Failed to create {:?}", &self.outdir))?;
        }
        let outdir = self
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
            &outdir,
            &gt_args,
        )?;
        info!("Successfully genotyped reads");
        // todo: create a filterer that takes a record and returns if the vcf record passes
        // todo: load pandora - passed - variants into HashMap with gene/interval keys
        // todo: load panel VCF and check with intersecting pandora vcf records
        let mut panel_vcf = rust_htslib::bcf::Reader::from_path(&self.index_vcf_path())
            .context("Failed to open panel VCF")?;
        for (i, record_result) in panel_vcf.records().enumerate() {
            let panel_variant = record_result
                .context(format!("Failed to read record {} in panel VCF", i))?;
        }
        // todo: generate prediction from pandora map output

        Ok(())
    }
}

impl Predict {
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

    fn validate_index(&self) -> Result<(), PredictError> {
        let expected_paths = &[
            self.index_prg_path(),
            self.index_kmer_prgs_path(),
            self.index_vcf_path(),
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
}

#[cfg(test)]
mod tests {
    use std::fs::File;

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
            let _f = File::create(tmp_path.join("genes.fa")).unwrap();
            let _f = File::create(tmp_path.join("dr.prg.k15.w14.idx")).unwrap();
        }

        let actual = predictor.validate_index().unwrap_err();
        let expected = PredictError::InvalidIndex(tmp_path.join("dr.update_DS"));

        assert_eq!(actual, expected)
    }
}
