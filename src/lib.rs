use std::path::{Path, PathBuf};

use log::debug;
use thiserror::Error;

const MAKE_PRG_BIN: &str = "make_prg";

/// A collection of custom errors relating to the working with files for this package.
#[derive(Error, Debug)]
pub enum MissingDependencies {
    /// Indicates that the path provided is not executable
    #[error("{0} is not executable or in $PATH")]
    NotExecutable(String),
}

pub struct MakePrg {
    executable: String,
}

impl MakePrg {
    /// Creates a MakePrg object from a path or from default if None given
    pub fn from_path(path: &Option<PathBuf>) -> Result<MakePrg, MissingDependencies> {
        let default = dependency_dir().join(MAKE_PRG_BIN);
        let executable = from_path_or(&path, &default);
        match (path, executable) {
            (_, Some(exec)) => Ok(MakePrg { executable: exec }),
            (Some(p), None) => Err(MissingDependencies::NotExecutable(String::from(
                p.to_string_lossy(),
            ))),
            (None, None) => {
                Err(MissingDependencies::NotExecutable(MAKE_PRG_BIN.to_owned()))
            }
        }
    }
}

pub struct MultipleSeqAligner {
    executable: String,
}

/// Check if an (optional) path is executable, and return it as a String. If no path is given, test
/// if the given default (or its file name) is executable and return it as a String if it is.
fn from_path_or(path: &Option<PathBuf>, default: &PathBuf) -> Option<String> {
    match path {
        Some(p) => {
            let executable = String::from(p.to_string_lossy());
            if is_executable(&executable) {
                Some(executable)
            } else {
                None
            }
        }
        None => {
            let default_fname = default.file_name()?.to_string_lossy();
            debug!(
                "No {} executable given. Trying default locations...",
                default_fname
            );
            if is_executable(&default.to_string_lossy()) {
                debug!("Found make_prg at {}", &default.to_string_lossy());
                Some(String::from(default.to_string_lossy()))
            } else if is_executable(&*default_fname) {
                debug!("Found {} on PATH", default_fname);
                Some(String::from(default_fname))
            } else {
                None
            }
        }
    }
}

pub fn dependency_dir() -> PathBuf {
    std::env::current_exe()
        .unwrap()
        .parent()
        .unwrap()
        .join("../../src/ext")
}

/// Checks whether the program is executable. If it is, it returns the full path to the
/// executable file
pub fn is_executable(program: &str) -> bool {
    let cmd = format!("command -v {}", program);
    let result = std::process::Command::new("sh")
        .args(&["-c", &cmd])
        .output();
    match result {
        Ok(output) => output.status.success(),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn path_is_executable() {
        let program = "ls";
        assert!(is_executable(program))
    }

    #[test]
    fn path_is_not_executable() {
        let program = "foobar";
        assert!(!is_executable(program))
    }

    #[test]
    fn from_path_or_where_path_is_executable() {
        let path = Some(PathBuf::from("/bin/ls"));
        let default = PathBuf::from("DEFAULT");
        let actual = from_path_or(&path, &default).unwrap();
        assert_eq!(actual, String::from(path.unwrap().to_string_lossy()))
    }

    #[test]
    fn from_path_or_where_path_isnt_executable_and_neither_is_default() {
        let path = Some(PathBuf::from("/bin/XZY"));
        let default = PathBuf::from("DEFAULT");
        let actual = from_path_or(&path, &default);
        assert!(actual.is_none())
    }

    #[test]
    fn from_path_or_where_path_is_none_and_default_isnt_executable() {
        let path = None;
        let default = PathBuf::from("DEFAULT");
        let actual = from_path_or(&path, &default);
        assert!(actual.is_none())
    }

    #[test]
    fn from_path_or_where_path_is_none_and_default_is_executable() {
        let path = None;
        let default = PathBuf::from("/bin/ls");
        let actual = from_path_or(&path, &default).unwrap();
        let expected = String::from(default.to_str().unwrap());
        assert_eq!(actual, expected)
    }

    #[test]
    fn from_path_or_where_path_is_none_default_isnt_executable_but_default_filename_is()
    {
        let path = None;
        let default = PathBuf::from("/XZY/ls");
        let actual = from_path_or(&path, &default).unwrap();
        let expected = String::from("ls");
        assert_eq!(actual, expected)
    }
}
