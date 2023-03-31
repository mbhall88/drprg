use crate::Runner;
use anyhow::Result;
use clap::Parser;
use flate2::read::GzDecoder;
use reqwest::Url;
use std::io::Cursor;
use std::path::PathBuf;
use tar::Archive;
use crate::config::INDEX_CONFIG;

lazy_static! {
    // we unwrap here because $HOME not being set is extremely unlikely
    static ref DEFAULT_OUTDIR: PathBuf = PathBuf::from(format!("{}/.drprg/", std::env::var("HOME").unwrap()));
}

#[derive(Parser, Debug, Default)]
pub struct Index {
    /// Download a prebuilt index
    #[clap(short, long)]
    download: bool,
    /// The name/path of the index to interact with
    #[clap()]
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
}

impl Runner for Index {
    fn run(&mut self) -> Result<()> {
        let (species, version) = match self.name.split_once('@') {
            Some((s, v)) => (s.to_owned(), v.to_owned()),
            None => (self.name.to_owned(), "latest".to_string())
        };

        // Download the tar archive file
        let mut response = reqwest::blocking::get(url)?;
        let mut buf = vec![];
        response.copy_to(&mut buf)?;

        // Decompress the tar archive file in memory
        let cursor = Cursor::new(buf);
        let tar = GzDecoder::new(cursor);
        let mut archive = Archive::new(tar);
        archive.unpack(&self.outdir)?;

        Ok(())
    }
}
