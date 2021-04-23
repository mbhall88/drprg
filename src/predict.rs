use crate::Runner;
use anyhow::Result;
use log::debug;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
pub struct Predict {}

impl Runner for Predict {
    fn run(&self) -> Result<()> {
        debug!("Predicting...");
        Ok(())
    }
}
