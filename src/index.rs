use crate::index::DownloadError::GitTreeFailed;
use crate::Runner;
use anyhow::{Context, Result};
use clap::Parser;
use drprg::unwrap_or_continue;
use flate2::read::GzDecoder;
use log::{debug, info};
use prettytable::{Row, Table};
use regex::Regex;
use serde_derive::Deserialize;
use std::collections::BTreeMap;
use std::io::Cursor;
use std::path::{Path, PathBuf};
use tar::Archive;
use thiserror::Error;

lazy_static! {
    // we unwrap here because $HOME not being set is extremely unlikely
    pub static ref DEFAULT_OUTDIR: PathBuf = PathBuf::from(format!("{}/.drprg/", std::env::var("HOME").unwrap()));
    static ref SPECIES_REGEX: Regex =
        Regex::new(r"^species/(?P<species1>\w+)/(?P<species2>\w+)-(?P<version>\w+)\.tar\.gz$").unwrap();
    static ref NUCLEOTIDES: Vec<&'static [u8]> = vec![b"A", b"C", b"G", b"T"];
}

type GitTree = BTreeMap<String, BTreeMap<(String, String), String>>;

#[derive(Error, Debug)]
pub enum DownloadError {
    #[error("Failed to find version {version} for species {species}")]
    UnknownVersion { version: String, species: String },
    #[error("Failed to retrieve the Git tree from GitHub with message {0}")]
    GitTreeFailed(String),
}

#[derive(Parser, Debug, Default)]
pub struct Index {
    /// Download a prebuilt index
    #[clap(short, long, conflicts_with = "list")]
    download: bool,
    /// List all available (and downloaded) indices
    #[clap(short, long)]
    list: bool,
    /// The name/path of the index to download
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

#[derive(Deserialize, Debug)]
#[allow(dead_code)]
struct Commit {
    sha: String,
    url: String,
    tree: Vec<SubTree>,
    truncated: bool,
}

#[derive(Deserialize, Debug)]
#[allow(dead_code)]
struct SubTree {
    path: String,
    mode: String,
    #[serde(rename = "type")]
    object_type: String,
    sha: String,
    size: Option<u64>,
    url: String,
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
        } else if self.list {
            list_indices(&self.outdir)?;
        }

        Ok(())
    }
}

fn load_available_indices() -> Result<GitTree> {
    let mut headers = reqwest::header::HeaderMap::new();
    headers.insert(
        "X-GitHub-Api-Version",
        reqwest::header::HeaderValue::from_static("2022-11-28"),
    );
    headers.insert(
        "Accept",
        reqwest::header::HeaderValue::from_static("application/vnd.github+json"),
    );
    headers.insert(
        "User-Agent",
        reqwest::header::HeaderValue::from_static("drprg"),
    );

    let client = reqwest::blocking::Client::builder()
        .connection_verbose(false)
        .default_headers(headers)
        .build()?;
    let res = client.get("https://api.github.com/repos/mbhall88/drprg-index/git/trees/main?recursive=1").send()?;

    if !res.status().is_success() {
        return Err(GitTreeFailed(res.status().to_string()).into());
    }
    let content: Commit = res.json()?;

    let mut config: GitTree = BTreeMap::new();

    for subtree in content.tree {
        let caps = SPECIES_REGEX.captures(&subtree.path);
        if let Some(captures) = caps {
            let species1 =
                unwrap_or_continue!(captures.name("species1").ok_or("no species"))
                    .as_str()
                    .to_string();
            let species2 =
                unwrap_or_continue!(captures.name("species2").ok_or("no species"))
                    .as_str()
                    .to_string();
            let version =
                unwrap_or_continue!(captures.name("version").ok_or("no version"))
                    .as_str()
                    .to_string();
            let url = format!("https://github.com/mbhall88/drprg-index/raw/main/species/{species1}/{species2}-{version}.tar.gz");

            let entry = config.entry(species1).or_default();
            entry.insert((version, species2), url);
        }
    }

    Ok(config)
}

fn download_indices(
    species: &str,
    version: &str,
    outdir: &Path,
    force: bool,
) -> Result<()> {
    let index_config =
        load_available_indices().context("Failed to load the available indices")?;
    for (spec, spec_conf) in index_config {
        if spec == species || species == "all" {
            let ((ver, spec2), url) = if version == "latest" {
                spec_conf.last_key_value()
            } else {
                spec_conf.get_key_value(&(version.to_string(), species.to_string()))
            }
            .ok_or_else(|| DownloadError::UnknownVersion {
                version: version.to_string(),
                species: spec.to_string(),
            })?;
            let outpath = outdir.join(&spec).join(format!("{spec2}-{ver}"));

            if outpath.exists() {
                if force {
                    debug!("{outpath:?} already exists. Removing it...");
                    std::fs::remove_dir_all(&outpath)?;
                } else {
                    info!("{spec} index version {ver} already downloaded. Skipping...");
                    continue;
                }
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

fn list_indices(outdir: &Path) -> Result<()> {
    let index_config =
        load_available_indices().context("Failed to load the available indices")?;
    let mut avail_tbl = Table::new();

    let is_verbose = log::max_level() == log::Level::Debug;
    let mut header = vec!["Name", "Species", "Version", "Downloaded"];
    if is_verbose {
        header.push("URL");
    }

    avail_tbl.add_row(Row::from(header));

    for (species, spec_conf) in index_config {
        for ((version, species2), url) in spec_conf {
            let mut row = vec![
                format!("{species}@{version}"),
                species.to_string(),
                version.to_string(),
            ];

            let outpath = outdir.join(&species).join(format!("{species2}-{version}"));
            if outpath.exists() {
                row.push("Y".to_string());
            } else {
                row.push("N".to_string());
            }

            if is_verbose {
                row.push(url.to_string());
            }
            avail_tbl.add_row(Row::from(row));
        }
    }

    avail_tbl.printstd();
    Ok(())
}
