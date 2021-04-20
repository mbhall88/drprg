use crate::Runner;
use anyhow::Result;
use drprg::{dependency_dir, MakePrg, MissingDependencies};
use log::{debug, info};
use std::path::PathBuf;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
pub struct Build {
    /// Path to pandora executable. Will try in src/ext or $PATH if not given
    #[structopt(short = "p", long = "pandora", parse(from_os_str), global = true)]
    pub pandora_exec: Option<PathBuf>,
    /// Path to make_prg executable. Will try in src/ext or $PATH if not given
    #[structopt(short = "m", long = "makeprg", parse(from_os_str))]
    pub makeprg_exec: Option<PathBuf>,
}

impl Runner for Build {
    fn run(&self) -> Result<()> {
        debug!("{:?}", &self);
        debug!("Building...");
        let makeprg = MakePrg::from_arg(&self.makeprg_exec)?;
        Ok(())
    }
}
