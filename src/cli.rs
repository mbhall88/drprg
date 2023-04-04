use clap::Parser;

use crate::builder::Build;
use crate::index::{Index, DEFAULT_OUTDIR};
use crate::predict::Predict;
use drprg::has_single_component;
use std::ffi::OsStr;
use std::path::PathBuf;

/// A utility function that allows the CLI to error if a path doesn't exist
pub fn check_path_exists<S: AsRef<OsStr> + ?Sized>(s: &S) -> Result<PathBuf, String> {
    let path = PathBuf::from(s);
    if path.exists() {
        Ok(path)
    } else {
        Err(format!("{path:?} does not exist"))
    }
}

/// A utility function to parse and validate the index provided to predict
pub fn validate_index(s: &str) -> Result<PathBuf, String> {
    let p = PathBuf::from(s);

    if p.exists() {
        Ok(p)
    } else if !has_single_component(s) {
        Err("Received an index which is path-like but does not exist".to_string())
    } else {
        let (species, version) = match s.split_once('@') {
            Some((sp, v)) => (sp.to_string(), v.to_string()),
            None => (s.to_owned(), "latest".to_string()),
        };

        let default_idx_dir = DEFAULT_OUTDIR.to_owned();
        let mut idx_path = default_idx_dir.join(&species);
        if !idx_path.exists() {
            return Err(format!(
                "No index for species {species} found in {default_idx_dir:?}"
            ));
        }
        match version.as_str() {
            "latest" => {
                let mut entries = std::fs::read_dir(&idx_path)
                    .map_err(|e| e.to_string())?
                    .map(|res| res.map(|e| e.path()))
                    .collect::<Result<Vec<_>, std::io::Error>>()
                    .map_err(|e| e.to_string())?;

                entries.sort_unstable();
                let mut d = entries.pop();

                while d.is_some() && !d.as_ref().unwrap().is_dir() {
                    d = entries.pop();
                }
                match d {
                    None => {
                        return Err(format!(
                        "No index versions found in {species} directory {idx_path:?}"
                    ))
                    }
                    Some(dir) => {
                        idx_path = dir;
                    }
                }
            }
            v => {
                let name = format!("{species}-{v}");
                idx_path.push(name);
                if !idx_path.exists() {
                    return Err(format!(
                        "Version {v} does not exist for species {species}"
                    ));
                }
            }
        }
        Ok(idx_path)
    }
}

/// Drug Resistance Prediction with Reference Graphs
#[derive(Parser, Debug)]
#[clap(author, version, about)]
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
    /// Build an index to predict resistance from
    Build(Build),
    /// Predict drug resistance
    Predict(Predict),
    /// Download and interact with indices
    Index(Index),
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
}
