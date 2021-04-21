use crate::Runner;
use anyhow::Result;
use drprg::{dependency_dir, MakePrg, MissingDependencies};
use log::{debug, info};
use std::path::PathBuf;
use structopt::StructOpt;
use bio::io::fasta;
use std::process::id;

#[derive(StructOpt, Debug)]
pub struct Build {
    /// Path to pandora executable. Will try in src/ext or $PATH if not given
    #[structopt(short = "p", long = "pandora", parse(from_os_str), global = true)]
    pub pandora_exec: Option<PathBuf>,
    /// Path to make_prg executable. Will try in src/ext or $PATH if not given
    #[structopt(short = "m", long = "makeprg", parse(from_os_str))]
    pub makeprg_exec: Option<PathBuf>,
    /// Annotation file that will be used to gather information about genes in panel
    #[structopt(short = "a", long = "gff", parse(from_os_str))]
    pub gff_file: PathBuf,
    /// Panel to build index from
    #[structopt(short = "i", long = "panel", parse(from_os_str))]
    pub panel_file: PathBuf,
    /// Reference genome in FASTA format (must be indexed with samtools faidx)
    #[structopt(short = "f", long = "fasta", parse(from_os_str))]
    pub reference_file: PathBuf,
    /// Number of bases of padding to add to start and end of each gene
    #[structopt(short = "P", long = "padding", default_value = "100")]
    pub padding: u32,
    /// Directory to place output
    #[structopt(short = "o", long = "outdir", default_value = ".")]
    pub outdir: PathBuf,
}

impl Runner for Build {
    fn run(&self) -> Result<()> {
        debug!("{:?}", &self);
        debug!("Building...");

        // get a set of all genes in the panel
        // load fasta index
        let mut faidx = fasta::IndexedReader::from_file(&self.reference_file)?;
        // open vcf ref file handle
        // iterate over gff and create hashmap of gene name to gff record
        // at the same time as above, extract the ref sequence for each gene

        // faidx.fetch_all("NC_000962.3")?;
        // let mut seq: Vec<u8> = vec![];
        // faidx.read(&mut seq);
        // println!("{:?}", seq.into_iter().map(|c| c as char).collect::<String>());


        // let makeprg = MakePrg::from_arg(&self.makeprg_exec)?;
        Ok(())
    }
}
