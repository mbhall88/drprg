use crate::Runner;
use anyhow::Result;
use log::debug;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
pub struct Build {}

impl Runner for Build {
    fn run(&self) -> Result<()> {
        debug!("Building...");
        Ok(())
    }
}
