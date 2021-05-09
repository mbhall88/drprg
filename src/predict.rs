use anyhow::{Context, Result};
use log::info;
use std::path::PathBuf;
use structopt::StructOpt;

use crate::cli::check_path_exists;
use crate::Runner;

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

        // todo: check all necessary index files exist
        Ok(())
    }
}
