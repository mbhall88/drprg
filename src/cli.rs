use structopt::StructOpt;

use crate::builder::Build;
use crate::predict::Predict;

// https://docs.rs/structopt/0.3.21/structopt/#subcommands
/// Drug Resistance Prediction with Reference Graphs
#[derive(StructOpt, Debug)]
pub struct Cli {
    /// Use verbose output
    #[structopt(short, long, global = true)]
    pub verbose: bool,
    #[structopt(subcommand)] // Note that we mark a field as a subcommand
    pub(crate) cmd: Command,
}

#[derive(StructOpt, Debug)]
pub enum Command {
    /// Build a drprg panel from a mykrobe-style panel
    Build(Build),
    /// Predict drug resistance
    Predict(Predict),
}
