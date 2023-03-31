use reqwest::Url;
use serde_derive::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

lazy_static! {
    pub static ref INDEX_CONFIG: HashMap<&'static str, HashMap<&'static str, Url>> = {
        let mut conf = HashMap::new();
        let species_versions = [("mtb", ["20230308"])];
        for (species, versions) in species_versions {
            let mut species_conf = HashMap::new();
            for ver in versions {
                let urlstr = format!("https://github.com/mbhall88/drprg-index/raw/main/species/{species}/{species}-{ver}.tar.gz");
                let url = Url::parse(&urlstr).unwrap();
                species_conf.insert(ver, url);
            }
            conf.insert(species, species_conf);
        }
        conf
    };
}

/// A collection of custom errors relating to the working with files for this package.
#[derive(Error, Debug)]
pub enum ConfigError {
    /// Generic I/O error for external dependencies
    #[error("File I/O failed in external dependency")]
    FileError { source: std::io::Error },
    /// Error trying to convert the TOML file to a string
    #[error(transparent)]
    DeserializeError(#[from] toml::de::Error),
}

#[derive(Deserialize, Serialize)]
pub struct Config {
    pub(crate) min_match_len: u32,
    pub(crate) max_nesting: u32,
    pub(crate) k: u32,
    pub(crate) w: u32,
    pub(crate) padding: u32,
    pub(crate) version: String,
}

impl Config {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, ConfigError> {
        let contents = std::fs::read_to_string(path)
            .map_err(|source| ConfigError::FileError { source })?;
        toml::from_str(&contents).map_err(ConfigError::DeserializeError)
    }
}
