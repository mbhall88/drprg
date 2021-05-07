use std::path::PathBuf;

use anyhow::Result;
use log::debug;
use structopt::StructOpt;

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
    #[structopt(short = "x", long, required = true, parse(from_os_str))]
    index: PathBuf,
    /// Sample reads to predict resistance from
    ///
    /// Both fasta and fastq are accepted, along with compressed or uncompressed.
    #[structopt(short, long, required = true, parse(from_os_str))]
    input: PathBuf,
    /// Directory to place output
    #[structopt(short, long, default_value = ".", parse(from_os_str))]
    outdir: PathBuf,
}

impl Runner for Predict {
    fn run(&self) -> Result<()> {
        debug!("Predicting...");
        Ok(())
    }
}
