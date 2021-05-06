use std::path::{Path, PathBuf};

use bstr::ByteSlice;
use log::{debug, error};
use std::ffi::OsStr;
use std::fs::File;
use std::process::Command;
use thiserror::Error;

const MAKE_PRG_BIN: &str = "make_prg";
const MAFFT_BIN: &str = "mafft/bin/mafft";
const PANDORA_BIN: &str = "pandora";

/// A collection of custom errors relating to the working with files for this package.
#[derive(Error, Debug)]
pub enum DependencyError {
    /// Indicates that the path provided is not executable
    #[error("{0} is not executable or in $PATH")]
    NotExecutable(String),
    /// Generic I/O error for external dependencies
    #[error("File I/O failed in external dependency")]
    FileError { source: std::io::Error },
    /// Represents all other cases of `std::io::Error`.
    #[error(transparent)]
    ProcessError(#[from] std::io::Error),
}

pub struct MakePrg {
    executable: String,
}

impl MakePrg {
    /// Creates a MakePrg object from a path or from default if None given
    pub fn from_path(path: &Option<PathBuf>) -> Result<MakePrg, DependencyError> {
        let default = dependency_dir().join(MAKE_PRG_BIN);
        let executable = from_path_or(&path, &default);
        match (path, executable) {
            (_, Some(exec)) => Ok(Self { executable: exec }),
            (Some(p), None) => Err(DependencyError::NotExecutable(String::from(
                p.to_string_lossy(),
            ))),
            (None, None) => {
                Err(DependencyError::NotExecutable(MAKE_PRG_BIN.to_owned()))
            }
        }
    }
    /// Run make_prg with the provided input, output and additional arguments
    pub fn run_with<I, S>(
        &self,
        input: &Path,
        output: &Path,
        args: I,
    ) -> Result<(), DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let dir = tempfile::tempdir().map_err(DependencyError::ProcessError)?;
        let prefix = "dr";
        let cmd_output = Command::new(&self.executable)
            .current_dir(&dir)
            .arg("from_msa")
            .args(args)
            .args(&["-o", prefix, "-i"])
            .arg(input)
            .output()
            .map_err(DependencyError::ProcessError)?;

        if !cmd_output.status.success() {
            error!(
                "Failed to run make_prg with sterr:\n{}",
                cmd_output.stderr.to_str_lossy()
            );
            Err(DependencyError::ProcessError(
                std::io::Error::from_raw_os_error(
                    cmd_output.status.code().unwrap_or(129),
                ),
            ))
        } else {
            debug!("make_prg successfully ran. Cleaning up temporary files...");
            let tmp_prefix = &dir.path().join(prefix);
            let tmp_prg = tmp_prefix.with_extension("prg.fa");
            let tmp_update_ds = tmp_prefix.with_extension("update_DS");

            // have to use fs::copy here as fs::rename fails inside a container as the tempdir is
            // not on the same "mount" as the local filesystem see more info at
            // https://doc.rust-lang.org/std/fs/fn.rename.html#platform-specific-behavior
            std::fs::copy(tmp_prg, output)
                .map_err(|source| DependencyError::FileError { source })?;
            std::fs::copy(tmp_update_ds, output.with_extension("update_DS"))
                .map_err(|source| DependencyError::FileError { source })?;

            Ok(())
        }
    }
}

pub struct Pandora {
    executable: String,
}

impl Pandora {
    pub fn from_path(path: &Option<PathBuf>) -> Result<Pandora, DependencyError> {
        let default = dependency_dir().join(PANDORA_BIN);
        let executable = from_path_or(&path, &default);
        match (path, executable) {
            (_, Some(exec)) => Ok(Self { executable: exec }),
            (Some(p), None) => Err(DependencyError::NotExecutable(String::from(
                p.to_string_lossy(),
            ))),
            (None, None) => Err(DependencyError::NotExecutable(
                default.file_name().unwrap().to_string_lossy().to_string(),
            )),
        }
    }

    /// Run pandora index with the provided input and arguments
    pub fn index_with<I, S>(&self, input: &Path, args: I) -> Result<(), DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let cmd_output = Command::new(&self.executable)
            .arg("index")
            .args(args)
            .arg(input)
            .output()
            .map_err(DependencyError::ProcessError)?;

        if !cmd_output.status.success() {
            error!(
                "Failed to run pandora index with sterr:\n{}",
                cmd_output.stderr.to_str_lossy()
            );
            Err(DependencyError::ProcessError(
                std::io::Error::from_raw_os_error(
                    cmd_output.status.code().unwrap_or(129),
                ),
            ))
        } else {
            Ok(())
        }
    }
}

pub struct MultipleSeqAligner {
    executable: String,
}

impl MultipleSeqAligner {
    pub fn from_path(
        path: &Option<PathBuf>,
    ) -> Result<MultipleSeqAligner, DependencyError> {
        let default = dependency_dir().join(MAFFT_BIN);
        let executable = from_path_or(&path, &default);
        match (path, executable) {
            (_, Some(exec)) => Ok(Self { executable: exec }),
            (Some(p), None) => Err(DependencyError::NotExecutable(String::from(
                p.to_string_lossy(),
            ))),
            (None, None) => Err(DependencyError::NotExecutable(
                default.file_name().unwrap().to_string_lossy().to_string(),
            )),
        }
    }
    /// Run the multiple sequence aligner with the provided input, output and arguments
    pub fn run_with<I, S>(
        &self,
        input: &Path,
        output: &Path,
        args: I,
    ) -> Result<(), DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let ostream = File::create(output)
            .map_err(|source| DependencyError::FileError { source })?;
        let result = Command::new(&self.executable)
            .args(args)
            .arg(input)
            .stdout(ostream)
            .output();
        match result {
            Ok(out) if out.status.success() => Ok(()),
            Ok(out) => {
                error!(
                    "Failed to run MAFFT with sterr:\n{}",
                    out.stderr.to_str_lossy()
                );
                Err(DependencyError::ProcessError(
                    std::io::Error::from_raw_os_error(out.status.code().unwrap_or(129)),
                ))
            }
            Err(e) => {
                error!("Failed to run MAFFT with sterr:\n");
                Err(DependencyError::ProcessError(e))
            }
        }
    }
}

/// Check if an (optional) path is executable, and return it as a String. If no path is given, test
/// if the given default (or its file name) is executable and return it as a String if it is.
fn from_path_or(path: &Option<PathBuf>, default: &Path) -> Option<String> {
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
                debug!("Found {} at {}", default_fname, &default.to_string_lossy());
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
    let result = Command::new("sh").args(&["-c", &cmd]).output();
    match result {
        Ok(output) => output.status.success(),
        _ => false,
    }
}

pub trait PathExt {
    fn add_extension(&self, extension: &OsStr) -> PathBuf;
}

impl PathExt for Path {
    fn add_extension(&self, extension: &OsStr) -> PathBuf {
        let mut s = self.as_os_str().to_os_string();
        s.push(extension);
        PathBuf::from(s)
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

    #[test]
    fn add_extension_empty_input() {
        let path = Path::new("foo.bar");

        let actual = path.add_extension("".as_ref());

        assert_eq!(actual, path)
    }

    #[test]
    fn add_extension_one_extension() {
        let path = Path::new("foo.bar");

        let actual = path.add_extension(".baz".as_ref());
        let expected = PathBuf::from("foo.bar.baz");

        assert_eq!(actual, expected)
    }

    #[test]
    fn add_extension_no_extension() {
        let path = Path::new("foo");

        let actual = path.add_extension(".baz".as_ref());
        let expected = PathBuf::from("foo.baz");

        assert_eq!(actual, expected)
    }
}
