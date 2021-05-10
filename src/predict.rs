use anyhow::{Context, Result};
use log::info;
use std::path::PathBuf;
use structopt::StructOpt;
use thiserror::Error;

use crate::cli::check_path_exists;
use crate::Runner;

/// A collection of custom errors relating to the predict component of this package
#[derive(Error, Debug, PartialEq)]
pub enum PredictError {
    /// The index is not valid
    #[error("Index is not valid due to missing file {0:?}")]
    InvalidIndex(PathBuf),
}

#[derive(StructOpt, Debug)]
pub struct Predict {
    /// Path to pandora executable. Will try in src/ext or $PATH if not given
    #[structopt(
        short = "p",
        long = "pandora",
        parse(from_os_str),
        hidden_short_help = true
    )]
    pandora_exec: Option<PathBuf>,
    /// Directory containing the index (produced by `drprg build`)
    #[structopt(short = "x", long, required = true, parse(try_from_os_str = check_path_exists))]
    index: PathBuf,
    /// Sample reads to predict resistance from
    ///
    /// Both fasta and fastq are accepted, along with compressed or uncompressed.
    #[structopt(short, long, required = true, parse(try_from_os_str = check_path_exists))]
    input: PathBuf,
    /// Directory to place output
    #[structopt(short, long, default_value = ".", parse(from_os_str))]
    outdir: PathBuf,
}

impl Runner for Predict {
    fn run(&self) -> Result<()> {
        if !self.outdir.exists() {
            info!("Outdir doesn't exist...creating...");
            std::fs::create_dir(&self.outdir)
                .context(format!("Failed to create {:?}", &self.outdir))?;
        }
        let _outdir = self
            .outdir
            .canonicalize()
            .context("Failed to canonicalize outdir")?;

        self.validate_index()?;
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
    use super::*;
    use std::fs::File;

    #[test]
    fn index_prg_path() {
        let predictor = Predict {
            pandora_exec: None,
            index: PathBuf::from("foo"),
            input: Default::default(),
            outdir: Default::default(),
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
        };

        let actual = predictor.index_prg_index_path();
        let expected = PathBuf::from("foo/dr.prg.k15.w14.idx");

        assert_eq!(actual, expected)
    }

    #[test]
    fn index_prg_update_path() {
        let predictor = Predict {
            pandora_exec: None,
            index: PathBuf::from("foo"),
            input: Default::default(),
            outdir: Default::default(),
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
            outdir: Default::default(),
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
            outdir: Default::default(),
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
            input: Default::default(),
            outdir: Default::default(),
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
