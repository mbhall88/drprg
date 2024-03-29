#[macro_use]
extern crate lazy_static;

use anyhow::Result;
use clap::Parser;
use env_logger::Builder;
use log::LevelFilter;
use log::{debug, info};

pub use crate::cli::Cli;
use crate::cli::Command;

pub mod builder;
mod cli;
pub mod config;
mod consequence;
mod expert;
mod index;
mod minor;
mod panel;
pub mod predict;
pub mod report;

pub trait Runner {
    fn run(&mut self) -> Result<()>;
}

fn main() -> Result<()> {
    let args = Cli::parse();
    // setup logging
    let log_lvl = if args.verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };
    let mut log_builder = Builder::new();
    log_builder
        .filter(None, log_lvl)
        .filter_module("reqwest", LevelFilter::Off)
        .format_module_path(false)
        .format_target(false)
        .init();
    debug!("{:?}", args);

    // set the global default number of threads for rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads as usize)
        .build_global()
        .unwrap();

    let mut subcmd: Box<dyn Runner> = match args.cmd {
        Command::Predict(cmd) => Box::new(cmd),
        Command::Build(cmd) => Box::new(cmd),
        Command::Index(cmd) => Box::new(cmd),
    };

    subcmd.run()?;

    info!("Done!");
    Ok(())
}
