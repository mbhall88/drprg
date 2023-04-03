use crate::config::INDEX_CONFIG;
use crate::Runner;
use anyhow::{Context, Result};
use clap::Parser;
use flate2::read::GzDecoder;
use log::{debug, info};
use std::io::Cursor;
use std::path::{Path, PathBuf};
use tar::Archive;
use thiserror::Error;

lazy_static! {
    // we unwrap here because $HOME not being set is extremely unlikely
    static ref DEFAULT_OUTDIR: PathBuf = PathBuf::from(format!("{}/.drprg/", std::env::var("HOME").unwrap()));
}

#[derive(Error, Debug)]
pub enum DownloadError {
    #[error("Failed to find version {version} for species {species}")]
    UnknownVersion { version: String, species: String },
}

#[derive(Parser, Debug, Default)]
pub struct Index {
    /// Download a prebuilt index
    #[clap(short, long, conflicts_with = "list")]
    download: bool,
    /// List all available (and downloaded) indices and version
    #[clap(short, long)]
    list: bool,
    /// The name/path of the index to interact with
    #[clap(default_value = "all")]
    name: String,
    /// Index directory
    ///
    /// Use this if your indices are not in a default location, or you want to download
    /// them to a non-default location
    #[clap(
    short,
    long,
    default_value_os_t = DEFAULT_OUTDIR.to_path_buf(),
    value_parser,
    value_name = "DIR"
    )]
    outdir: PathBuf,
    /// Overwrite any existing indices
    #[clap(short = 'F', long)]
    force: bool,
}

impl Runner for Index {
    fn run(&mut self) -> Result<()> {
        let (species, version) = match self.name.split_once('@') {
            Some((s, v)) => (s.to_owned(), v.to_owned()),
            None => (self.name.to_owned(), "latest".to_string()),
        };

        if self.download {
            download_indices(&species, &version, &self.outdir, self.force).context(
                format!("Failed to download {species} species version {version}"),
            )?;
        }

        Ok(())
    }
}

fn download_indices(
    species: &str,
    version: &str,
    outdir: &Path,
    force: bool,
) -> Result<()> {
    for (spec, spec_conf) in &*INDEX_CONFIG {
        if spec == &species || species == "all" {
            let (ver, url) = if version == "latest" {
                spec_conf.last_key_value()
            } else {
                spec_conf.get_key_value(version)
            }
            .ok_or_else(|| DownloadError::UnknownVersion {
                version: version.to_string(),
                species: spec.to_string(),
            })?;
            let outpath = outdir.join(spec).join(format!("{spec}-{ver}"));

            if outpath.exists() && force {
                debug!("{outpath:?} already exists. Removing it...");
                std::fs::remove_dir_all(&outpath)?;
            } else {
                info!("{spec} index version {ver} already downloaded. Skipping...");
                continue;
            }

            info!("Downloading {spec} index version {ver} to {outpath:?}...");
            download_from_url(url, outpath.parent().unwrap())?;
            info!("Download complete");
        }
    }

    Ok(())
}

fn download_from_url(url: &str, dest: &Path) -> Result<()> {
    // Download the tar archive file
    let mut response = reqwest::blocking::get(url)?;
    let mut buf = vec![];
    response.copy_to(&mut buf)?;

    // Decompress the tar archive file in memory
    let cursor = Cursor::new(buf);
    let tar = GzDecoder::new(cursor);
    let mut archive = Archive::new(tar);
    archive.unpack(dest)?;
    Ok(())
}
