use clap::Parser;

use crate::builder::Build;
use crate::predict::Predict;
use std::ffi::OsStr;
use std::path::PathBuf;
use drprg::PathExt;

/// A utility function that allows the CLI to error if a path doesn't exist
pub fn check_path_exists<S: AsRef<OsStr> + ?Sized>(s: &S) -> Result<PathBuf, String> {
    let path = PathBuf::from(s);
    if path.exists() {
        Ok(path)
    } else {
        Err(format!("{:?} does not exist", path))
    }
}

/// A utility function that allows the CLI to error if a VCF and its index doesn't exist
pub fn check_vcf_and_index_exists<S: AsRef<OsStr> + ?Sized>(s: &S) -> Result<PathBuf, String> {
    let path = PathBuf::from(s);
    if !&path.exists() {
        Err(format!("{:?} does not exist", &path))
    } else {
        for ext in &[".csi", ".tbi"] {
            let p = &path.add_extension(ext.as_ref());
            if p.exists() {
                return Ok(path);
            }
        }
        Err(format!("{:?} does not have an index", path))
    }
}

/// Drug Resistance Prediction with Reference Graphs
#[derive(Parser, Debug)]
pub struct Cli {
    /// Use verbose output
    #[clap(short, long, global = true)]
    pub verbose: bool,
    /// Maximum number of threads to use
    ///
    /// Use 0 to select the number automatically
    #[clap(short, long, global = true, default_value = "1", value_name = "INT")]
    pub threads: u8,
    #[clap(subcommand)] // Note that we mark a field as a subcommand
    pub(crate) cmd: Command,
}

#[derive(Parser, Debug)]
pub enum Command {
    /// Build a drprg panel from a mykrobe-style panel
    Build(Build),
    /// Predict drug resistance
    Predict(Predict),
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_path_exists_it_doesnt() {
        let result = check_path_exists(OsStr::new("fake.path"));
        assert!(result.is_err())
    }

    #[test]
    fn check_path_it_does() {
        let actual = check_path_exists(OsStr::new("Cargo.toml")).unwrap();
        let expected = PathBuf::from("Cargo.toml");
        assert_eq!(actual, expected)
    }

    #[test]
    fn check_vcf_and_index_exists_it_doesnt() {
        let actual = check_vcf_and_index_exists(OsStr::new("tests/cases/predict/in.vcf"));
        assert!(actual.is_err())
    }

    #[test]
    fn check_vcf_and_index_exists_it_does() {
        let path = OsStr::new("tests/cases/build/input.bcf");
        let actual = check_vcf_and_index_exists(path).unwrap();
        let expected = PathBuf::from(path);
        assert_eq!(actual, expected)
    }
}
