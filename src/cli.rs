use crate::build::build::Build;
use crate::predict::predict::Predict;
use structopt::StructOpt;

// https://docs.rs/structopt/0.3.21/structopt/#subcommands
/// Drug Resistance Prediction with Reference Graphs
#[derive(StructOpt, Debug)]
pub struct Cli {
    /// Use verbose output
    #[structopt(short, global = true)]
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
