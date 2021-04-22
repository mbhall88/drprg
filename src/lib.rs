use std::path::PathBuf;

use log::debug;
use thiserror::Error;

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
    /// Creates a MakePrg object from a CLI arg or from default if None given
    pub fn from_arg(arg: &Option<PathBuf>) -> Result<MakePrg, MissingDependencies> {
        match arg {
            Some(p) => {
                let executable = String::from(p.to_string_lossy());
                if is_executable(&executable) {
                    Ok(MakePrg { executable })
                } else {
                    Err(MissingDependencies::NotExecutable(executable))
                }
            }
            None => {
                debug!("No make_prg executable given. Trying default locations...");
                let executable =
                    String::from(dependency_dir().join("make_prg").to_string_lossy());
                if is_executable(&executable) {
                    debug!("Found make_prg at {}", &executable);
                    Ok(MakePrg { executable })
                } else if is_executable("make_prg") {
                    debug!("Found make_prg on PATH");
                    Ok(MakePrg {
                        executable: String::from("make_prg"),
                    })
                } else {
                    Err(MissingDependencies::NotExecutable(executable))
                }
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
}
