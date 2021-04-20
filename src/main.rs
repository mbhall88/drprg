use anyhow::Result;
use env_logger::Builder;
use log::LevelFilter;
use log::{debug, info};
use structopt::StructOpt;

pub use crate::cli::Cli;
use crate::cli::Command;

mod build;
mod cli;
mod predict;

pub trait Runner {
    fn run(&self) -> Result<()>;
}

fn main() -> Result<()> {
    let args: Cli = Cli::from_args();
    // setup logging
    let log_lvl = if args.verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };
    let mut log_builder = Builder::new();
    log_builder
        .filter(None, log_lvl)
        .format_module_path(false)
        .init();
    debug!("{:?}", args);

    let subcmd: Box<dyn Runner> = match args.cmd {
        Command::Predict(cmd) => Box::new(cmd),
        Command::Build(cmd) => Box::new(cmd),
    };

    subcmd.run()?;

    info!("Done!");
    Ok(())
}
