use std::ffi::OsStr;
use std::fs::File;
use std::ops::Range;
use std::os::unix::ffi::OsStringExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use bstr::ByteSlice;
use log::{debug, error};
use rust_htslib::bcf;
use rust_htslib::bcf::header::Id;
use rust_htslib::bcf::record::GenotypeAllele;
use thiserror::Error;

use crate::filter::Tags;
use crate::interval::IntervalOp;
use noodles::fasta::record::Definition;
use noodles::fasta::repository::adapters::IndexedReader;
use noodles::fasta::repository::Adapter;
use noodles::fasta::{fai, Record};
use noodles::{fasta, gff};
use regex::Regex;
use std::cmp::{max, min};
use std::fs;
use std::io::{BufReader, BufWriter, ErrorKind};

pub mod filter;
pub mod interval;

const MAKE_PRG_BIN: &str = "make_prg";
const MAFFT_BIN: &str = "mafft/bin/mafft";
const PANDORA_BIN: &str = "pandora";
const BCFTOOLS_BIN: &str = "bcftools";
const MTB_GENOME_SIZE: u32 = 4411532;

#[macro_export]
/// A macro that will unwrap and return the value of a Result if Ok or cause a loop to continue if
/// the Result is an Err
///
/// # Examples
/// ```rust
/// # use drprg::unwrap_or_continue;
/// let v = &[Ok(1), Ok(2), Err(3), Ok(4)];
/// let mut sum = 0;
/// let mut i = 0;
/// for x in v {
/// i += 1;
/// sum += unwrap_or_continue!(x);
/// }
/// assert_eq!(i, 4);
/// assert_eq!(sum, 7)
/// ```
macro_rules! unwrap_or_continue {
    ( $result:expr ) => {
        match $result {
            Ok(val) => val,
            Err(_) => continue,
        }
    };
}

/// A collection of custom errors relating to the working with files for this package.
#[derive(Error, Debug)]
pub enum DependencyError {
    /// Indicates that the path provided is not executable
    #[error("{0} is not executable or in $PATH")]
    NotExecutable(String),
    /// Generic I/O error for external dependencies
    #[error("File I/O failed in external dependency")]
    FileError { source: std::io::Error },
    /// Represents all other cases of `std::io::Error`.
    #[error(transparent)]
    ProcessError(#[from] std::io::Error),
    /// An error associated with indexing VCFs or fastas with htslib
    #[error("Failed to index VCF with htslib: {0}")]
    HtslibIndexError(String),
    /// Missing an expected output file after running a dependency
    #[error("Missing expected output file {0}")]
    MissingExpectedOutput(String),
    /// Denovo paths error
    #[error("Failed to extract information about novel variants: {0}")]
    NovelVariantParsingError(String),
    /// Issue parsing fasta records
    #[error("Failed to parse fasta record in {0:?}")]
    FastaParserError(PathBuf),
}

pub struct Bcftools {
    executable: String,
}

impl Bcftools {
    pub fn from_path(path: &Option<PathBuf>) -> Result<Self, DependencyError> {
        let default = dependency_dir().join(BCFTOOLS_BIN);
        let executable = from_path_or(path, &default);
        match (path, executable) {
            (_, Some(exec)) => Ok(Self { executable: exec }),
            (Some(p), None) => Err(DependencyError::NotExecutable(String::from(
                p.to_string_lossy(),
            ))),
            (None, None) => {
                Err(DependencyError::NotExecutable(BCFTOOLS_BIN.to_owned()))
            }
        }
    }

    pub fn sort(&self, input: &Path, output: &Path) -> Result<(), DependencyError> {
        let mut binding = Command::new(&self.executable);
        let cmd = binding
            .arg("sort")
            .args(["-O", "b", "-o"])
            .arg(output)
            .arg(input);
        let cmd_output = cmd.output().map_err(DependencyError::ProcessError)?;
        let cmd_args = cmd
            .get_args()
            .collect::<Vec<&OsStr>>()
            .join(OsStr::new(" "));
        debug!(
            "Ran: {} {}",
            cmd.get_program().to_string_lossy(),
            cmd_args.to_string_lossy()
        );

        if !cmd_output.status.success() {
            let exit_code = cmd_output.status.to_string();
            Err(DependencyError::ProcessError(std::io::Error::new(
                ErrorKind::Other,
                format!(
                    "{} with stderr: {}",
                    exit_code,
                    cmd_output.stderr.to_str_lossy()
                ),
            )))
        } else {
            Ok(())
        }
    }

    pub fn consensus<I, S>(
        &self,
        input: &Path,
        output: &Path,
        args: I,
    ) -> Result<(), DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let mut binding = Command::new(&self.executable);
        let cmd = binding
            .arg("consensus")
            .args(args)
            .arg("-o")
            .arg(output)
            .arg(input);
        let cmd_output = cmd.output().map_err(DependencyError::ProcessError)?;
        let cmd_args = cmd
            .get_args()
            .collect::<Vec<&OsStr>>()
            .join(OsStr::new(" "));
        debug!(
            "Ran: {} {}",
            cmd.get_program().to_string_lossy(),
            cmd_args.to_string_lossy()
        );

        if !cmd_output.status.success() {
            let exit_code = cmd_output.status.to_string();
            Err(DependencyError::ProcessError(std::io::Error::new(
                ErrorKind::Other,
                format!(
                    "{} with stderr: {}",
                    exit_code,
                    cmd_output.stderr.to_str_lossy()
                ),
            )))
        } else {
            Ok(())
        }
    }
}

pub struct MakePrg {
    executable: String,
}

impl MakePrg {
    /// Creates a MakePrg object from a path or from default if None given
    pub fn from_path(path: &Option<PathBuf>) -> Result<MakePrg, DependencyError> {
        let default = dependency_dir().join(MAKE_PRG_BIN);
        let executable = from_path_or(path, &default);
        match (path, executable) {
            (_, Some(exec)) => Ok(Self { executable: exec }),
            (Some(p), None) => Err(DependencyError::NotExecutable(String::from(
                p.to_string_lossy(),
            ))),
            (None, None) => {
                Err(DependencyError::NotExecutable(MAKE_PRG_BIN.to_owned()))
            }
        }
    }
    /// Run make_prg with the provided input, output and additional arguments
    pub fn from_msas_with<I, S>(
        &self,
        input: &Path,
        output_prg: &Path,
        output_update_prgs: &Path,
        args: I,
    ) -> Result<(), DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let dir = tempfile::tempdir().map_err(DependencyError::ProcessError)?;
        let prefix = "dr";
        let outdir = output_prg.parent().unwrap_or_else(|| Path::new("."));
        let logstream = File::create(outdir.join("makeprg.log"))
            .map_err(|source| DependencyError::FileError { source })?;
        let mut binding = Command::new(&self.executable);
        let cmd = binding
            .current_dir(&dir)
            .arg("from_msa")
            .args(args)
            .args(["-v", "-o", prefix, "-i"])
            .arg(input)
            .stdout(Stdio::null())
            .stderr(logstream);
        let cmd_result = cmd.output();
        let cmd_args = cmd
            .get_args()
            .collect::<Vec<&OsStr>>()
            .join(OsStr::new(" "));
        debug!(
            "Ran: {} {}",
            cmd.get_program().to_string_lossy(),
            cmd_args.to_string_lossy()
        );

        match cmd_result {
            Ok(cmd_output) if !cmd_output.status.success() => {
                let _ = fs::remove_dir_all(dir.path());
                error!(
                    "Failed to run make_prg with sterr:\n{}",
                    cmd_output.stderr.to_str_lossy()
                );
                Err(DependencyError::ProcessError(
                    std::io::Error::from_raw_os_error(
                        cmd_output.status.code().unwrap_or(129),
                    ),
                ))
            }
            Ok(_) => {
                debug!("make_prg successfully ran. Cleaning up temporary files...");
                let tmp_prefix = &dir.path().join(prefix);
                let tmp_prg = tmp_prefix.with_extension("prg.fa");
                let tmp_update_prgs = &dir.path().join(format!("{}_prgs", prefix));

                // have to use fs::copy here as fs::rename fails inside a container as the tempdir is
                // not on the same "mount" as the local filesystem see more info at
                // https://doc.rust-lang.org/std/fs/fn.rename.html#platform-specific-behavior
                fs::copy(tmp_prg, output_prg)
                    .map_err(|source| DependencyError::FileError { source })?;
                let copyopts = fs_extra::dir::CopyOptions {
                    overwrite: true,
                    skip_exist: false,
                    buffer_size: 64000,
                    copy_inside: true,
                    content_only: false,
                    depth: 0,
                };
                fs_extra::dir::move_dir(tmp_update_prgs, output_update_prgs, &copyopts)
                    .map_err(|source| DependencyError::FileError {
                        source: std::io::Error::new(
                            std::io::ErrorKind::Other,
                            source.to_string(),
                        ),
                    })?;

                let _ = fs::remove_dir_all(dir.path());
                Ok(())
            }
            Err(err) => {
                let _ = fs::remove_dir_all(dir.path());
                error!("make_prg failed to run with error: {}", err.to_string());
                Err(DependencyError::ProcessError(err))
            }
        }
    }

    /// Update a PRG with make_prg
    pub fn update<I, S>(
        &self,
        index_msas_dir: &Path,
        index_prg: &Path,
        denovo_paths: &Path,
        outdir: &Path,
        args: I,
        mafft: String,
    ) -> Result<PathBuf, DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let dir = outdir.join("tmp_msas");
        if !dir.exists() {
            fs::create_dir(&dir)?;
        }
        let genes_to_update = Pandora::list_prgs_with_novel_variants(denovo_paths)?;
        debug!(
            "{} genes to update with novel variants: {:?}",
            genes_to_update.len(),
            genes_to_update
        );
        let update_msa_dir = outdir.join("update_msas");
        if !update_msa_dir.exists() {
            fs::create_dir(&update_msa_dir)?;
        }
        let aligner = MultipleSeqAligner::from_path(&Some(PathBuf::from(&mafft)))?;
        let denovo_sequences =
            denovo_paths.parent().unwrap().join("denovo_sequences.fa");
        if !denovo_sequences.exists() {
            return Err(DependencyError::MissingExpectedOutput(
                denovo_sequences.to_string_lossy().to_string(),
            ));
        }

        let denovo_sequences_faidx = index_fasta(&denovo_sequences)?;
        if !denovo_sequences_faidx.exists() {
            return Err(DependencyError::MissingExpectedOutput(
                denovo_sequences_faidx.to_string_lossy().to_string(),
            ));
        }

        let fa_reader = File::open(&denovo_sequences)
            .map(BufReader::new)
            .map(fasta::Reader::new)?;
        let index = File::open(&denovo_sequences_faidx)
            .map(BufReader::new)
            .map(fai::Reader::new)?
            .read_index()?;
        let mut faidx = IndexedReader::new(fa_reader, index);

        for gene in &genes_to_update {
            debug!("Updating MSA for {}", gene);
            let existing_msa = index_msas_dir.join(format!("{}.fa", gene));
            let new_sequence = dir.join(format!("{}.consensus.fa", gene));
            {
                let record = match faidx.get(gene) {
                    Some(Ok(r)) => r,
                    _ => return Err(DependencyError::HtslibIndexError(format!("Could not extract {} from pandora discover consensus sequences", gene)))
                };
                let mut fa_writer =
                    File::create(&new_sequence).map(BufWriter::new).map(|f| {
                        fasta::Writer::builder(f)
                            .set_line_base_count(usize::MAX)
                            .build()
                    })?;
                let def = Definition::new(format!("{}_updated", gene), None);
                let new_record = Record::new(def, record.sequence().to_owned());
                fa_writer.write_record(&new_record)?;
            }
            let updated_msa = update_msa_dir.join(format!("{}.fa", gene));
            aligner.run_with(
                &existing_msa,
                &updated_msa,
                [
                    "--auto",
                    "--thread",
                    "-1",
                    "--quiet",
                    "--add",
                    new_sequence.to_str().unwrap(),
                ],
            )?;
        }

        let update_prgs_dir = outdir.join("update_prgs");
        if !update_prgs_dir.exists() {
            fs::create_dir(&update_prgs_dir)?;
        }
        let prefix = "updated";
        let logstream = File::create(update_prgs_dir.join("update_prgs.log"))
            .map_err(|source| DependencyError::FileError { source })?;

        let fixed_args = vec!["-v", "-o", prefix, "-i"];

        let mut binding = Command::new(&self.executable);
        let cmd = binding
            .current_dir(&update_prgs_dir)
            .arg("from_msa")
            .args(args)
            .args(&fixed_args)
            .arg(update_msa_dir)
            .stdout(Stdio::null())
            .stderr(logstream);
        let cmd_result = cmd.output();
        let cmd_args = cmd
            .get_args()
            .collect::<Vec<&OsStr>>()
            .join(OsStr::new(" "));
        debug!(
            "Ran: {} {}",
            cmd.get_program().to_string_lossy(),
            cmd_args.to_string_lossy()
        );

        match cmd_result {
            Ok(cmd_output) if !cmd_output.status.success() => {
                error!("Failed to run make_prg update. Check update_prgs.log",);
                Err(DependencyError::ProcessError(
                    std::io::Error::from_raw_os_error(
                        cmd_output.status.code().unwrap_or(129),
                    ),
                ))
            }
            Ok(_) => {
                debug!("make_prg update successfully ran");
                let updated_prgs =
                    update_prgs_dir.join(prefix).with_extension("prg.fa");
                {
                    let mut prg_writer = File::options()
                        .append(true)
                        .open(&updated_prgs)
                        .map(BufWriter::new)
                        .map(|f| {
                            fasta::Writer::builder(f)
                                .set_line_base_count(usize::MAX)
                                .build()
                        })?;
                    let mut prg_reader = File::open(index_prg)
                        .map(BufReader::new)
                        .map(fasta::Reader::new)?;

                    for res in prg_reader.records() {
                        let record = res.map_err(|_| {
                            DependencyError::FastaParserError(index_prg.to_path_buf())
                        })?;
                        let name = record.name().to_string();
                        if !genes_to_update.contains(&name) {
                            prg_writer.write_record(&record)?;
                        }
                    }
                }
                let output_prg = outdir.join("updated.dr.prg");
                fs::copy(&updated_prgs, &output_prg)
                    .map_err(|source| DependencyError::FileError { source })?;
                let _ = fs::remove_dir_all(&dir);
                Ok(output_prg)
            }
            Err(err) => {
                error!(
                    "make_prg update failed to run with error: {}",
                    err.to_string()
                );
                Err(DependencyError::ProcessError(err))
            }
        }
    }
}

pub struct Pandora {
    executable: String,
}

impl Pandora {
    pub fn from_path(path: &Option<PathBuf>) -> Result<Pandora, DependencyError> {
        let default = dependency_dir().join(PANDORA_BIN);
        let executable = from_path_or(path, &default);
        match (path, executable) {
            (_, Some(exec)) => Ok(Self { executable: exec }),
            (Some(p), None) => Err(DependencyError::NotExecutable(String::from(
                p.to_string_lossy(),
            ))),
            (None, None) => Err(DependencyError::NotExecutable(
                default.file_name().unwrap().to_string_lossy().to_string(),
            )),
        }
    }

    /// Run pandora index with the provided input and arguments
    pub fn index_with<I, S>(&self, input: &Path, args: I) -> Result<(), DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let mut binding = Command::new(&self.executable);
        let cmd = binding.arg("index").args(args).arg(input);
        let cmd_output = cmd.output().map_err(DependencyError::ProcessError)?;
        let cmd_args = cmd
            .get_args()
            .collect::<Vec<&OsStr>>()
            .join(OsStr::new(" "));
        debug!(
            "Ran: {} {}",
            cmd.get_program().to_string_lossy(),
            cmd_args.to_string_lossy()
        );

        if !cmd_output.status.success() {
            error!(
                "Failed to run pandora index with sterr:\n{}",
                cmd_output.stderr.to_str_lossy()
            );
            Err(DependencyError::ProcessError(
                std::io::Error::from_raw_os_error(
                    cmd_output.status.code().unwrap_or(129),
                ),
            ))
        } else {
            Ok(())
        }
    }

    /// Run pandora discover
    pub fn discover_with<I, S>(
        &self,
        prg: &Path,
        query_idx: &Path,
        outdir: &Path,
        args: I,
    ) -> Result<PathBuf, DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        if !outdir.exists() {
            std::fs::create_dir(outdir)?;
        }
        let logstream = File::create(outdir.join("discover.log"))
            .map_err(|source| DependencyError::FileError { source })?;
        let fixed_args = &["discover", "-g", &MTB_GENOME_SIZE.to_string(), "-v", "-o"];
        let mut binding = Command::new(&self.executable);
        let cmd = binding
            .args(fixed_args)
            .arg(outdir)
            .args(args)
            .arg(prg)
            .arg(query_idx)
            .stdout(logstream)
            .stderr(Stdio::inherit());
        let cmd_output = cmd.output().map_err(DependencyError::ProcessError)?;
        let cmd_args = cmd
            .get_args()
            .collect::<Vec<&OsStr>>()
            .join(OsStr::new(" "));
        debug!(
            "Ran: {} {}",
            cmd.get_program().to_string_lossy(),
            cmd_args.to_string_lossy()
        );

        if !cmd_output.status.success() {
            error!(
                "Failed to run pandora discover with sterr:\n{}",
                cmd_output.stderr.to_str_lossy()
            );
            Err(DependencyError::ProcessError(
                std::io::Error::from_raw_os_error(
                    cmd_output.status.code().unwrap_or(129),
                ),
            ))
        } else {
            let denovo_paths = outdir.join("denovo_paths.txt");
            if denovo_paths.exists() {
                Ok(denovo_paths)
            } else {
                Err(DependencyError::MissingExpectedOutput(
                    denovo_paths.to_string_lossy().to_string(),
                ))
            }
        }
    }

    pub fn genotype_with<I, S>(
        &self,
        prg: &Path,
        vcf_ref: &Path,
        reads: &Path,
        outdir: &Path,
        args: I,
    ) -> Result<(), DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let errstream = File::create(outdir.join("pandora.log"))
            .map_err(|source| DependencyError::FileError { source })?;
        let fixed_args = &[
            "map",
            "--genotype",
            "-v",
            "-o",
            &outdir.to_string_lossy(),
            "-g",
            &MTB_GENOME_SIZE.to_string(),
            "--vcf-refs",
            &vcf_ref.to_string_lossy(),
        ];
        let mut binding = Command::new(&self.executable);
        let cmd = binding
            .args(fixed_args)
            .args(args)
            .arg(prg)
            .arg(reads)
            .stdout(errstream)
            .stderr(Stdio::inherit());
        let cmd_output = cmd.output().map_err(DependencyError::ProcessError)?;
        let cmd_args = cmd
            .get_args()
            .collect::<Vec<&OsStr>>()
            .join(OsStr::new(" "));
        debug!(
            "Ran: {} {}",
            cmd.get_program().to_string_lossy(),
            cmd_args.to_string_lossy()
        );

        if !cmd_output.status.success() {
            error!(
                "Failed to run pandora map with sterr:\n{}",
                cmd_output.stderr.to_str_lossy()
            );
            Err(DependencyError::ProcessError(
                std::io::Error::from_raw_os_error(
                    cmd_output.status.code().unwrap_or(129),
                ),
            ))
        } else {
            Ok(())
        }
    }

    pub fn vcf_filename() -> String {
        String::from("pandora_genotyped.vcf")
    }

    pub fn list_prgs_with_novel_variants(
        denovo_file: &Path,
    ) -> Result<Vec<String>, DependencyError> {
        let mut genes = vec![];

        let re = Regex::new(r"\n(?P<num>\d+) loci with denovo variants\n").unwrap();
        let contents = std::fs::read_to_string(denovo_file).map_err(|_| {
            DependencyError::NovelVariantParsingError(format!(
                "Unable to read {:?}",
                denovo_file
            ))
        })?;
        let captures = re.captures(&contents).ok_or_else(|| {
            DependencyError::NovelVariantParsingError(
                "Unable to find line describing the number of novel variants"
                    .to_string(),
            )
        })?;
        // we unwrap here as we know we have a positive integer if the capture exists
        let expected_num_genes: usize = captures
            .name("num")
            .ok_or_else(|| {
                DependencyError::NovelVariantParsingError(
                    "Unable to find line describing the number of novel variants"
                        .to_string(),
                )
            })?
            .as_str()
            .parse()
            .unwrap();

        let mut previous_line = String::new();

        for line in contents.lines() {
            if line.ends_with("nodes") {
                genes.push(previous_line)
            }

            previous_line = line.to_string();
        }

        if genes.len() != expected_num_genes {
            Err(DependencyError::NovelVariantParsingError(format!(
                "Expected {} genes with novel variants, but found {}",
                expected_num_genes,
                genes.len()
            )))
        } else {
            Ok(genes)
        }
    }
}

pub struct MultipleSeqAligner {
    pub executable: String,
}

impl MultipleSeqAligner {
    pub fn from_path(
        path: &Option<PathBuf>,
    ) -> Result<MultipleSeqAligner, DependencyError> {
        let default = dependency_dir().join(MAFFT_BIN);
        let executable = from_path_or(path, &default);
        match (path, executable) {
            (_, Some(exec)) => Ok(Self { executable: exec }),
            (Some(p), None) => Err(DependencyError::NotExecutable(String::from(
                p.to_string_lossy(),
            ))),
            (None, None) => Err(DependencyError::NotExecutable(
                default.file_name().unwrap().to_string_lossy().to_string(),
            )),
        }
    }
    /// Run the multiple sequence aligner with the provided input, output and arguments
    pub fn run_with<I, S>(
        &self,
        input: &Path,
        output: &Path,
        args: I,
    ) -> Result<(), DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let ostream = File::create(output)
            .map_err(|source| DependencyError::FileError { source })?;
        let mut binding = Command::new(&self.executable);
        let cmd = binding.args(args).arg(input).stdout(ostream);
        let cmd_result = cmd.output();
        let cmd_args = cmd
            .get_args()
            .collect::<Vec<&OsStr>>()
            .join(OsStr::new(" "));
        debug!(
            "Ran: {} {}",
            cmd.get_program().to_string_lossy(),
            cmd_args.to_string_lossy()
        );
        match cmd_result {
            Ok(out) if out.status.success() => Ok(()),
            Ok(out) => {
                error!(
                    "Failed to run MAFFT with sterr:\n{}",
                    out.stderr.to_str_lossy()
                );
                Err(DependencyError::ProcessError(
                    std::io::Error::from_raw_os_error(out.status.code().unwrap_or(129)),
                ))
            }
            Err(e) => {
                error!("Failed to run MAFFT with sterr:\n");
                Err(DependencyError::ProcessError(e))
            }
        }
    }
}

/// Check if an (optional) path is executable, and return it as a String. If no path is given, test
/// if the given default (or its file name) is executable and return it as a String if it is.
fn from_path_or(path: &Option<PathBuf>, default: &Path) -> Option<String> {
    match path {
        Some(p) => {
            let executable = String::from(p.to_string_lossy());
            is_executable(&executable)
        }
        None => {
            let default_fname = default.file_name()?.to_string_lossy();
            is_executable(&default.to_string_lossy())
                .or_else(|| is_executable(&default_fname))
        }
    }
}

pub fn index_vcf(path: &Path) -> Result<(), DependencyError> {
    let min_shift: i32 = 14; // recommended default in htslib docs - https://github.com/samtools/htslib/blob/818008a750eefb347bb3732dff9fb60afc367de6/htslib/vcf.h#L1236
    let fname = std::ffi::CString::new(path.to_path_buf().into_os_string().into_vec())
        .map_err(|_| {
            DependencyError::HtslibIndexError(
                "Failed to convert filtered VCF path into a CString".to_string(),
            )
        })?;
    unsafe {
        match rust_htslib::htslib::bcf_index_build(fname.as_ptr(), min_shift) {
            0 => Ok(()),
            -1 => Err(DependencyError::HtslibIndexError(
                "Indexing failed (htslib exit code -1)".to_string(),
            )),
            -2 => Err(DependencyError::HtslibIndexError(
                "Opening @fn failed (htslib exit code -2)".to_string(),
            )),
            -3 => Err(DependencyError::HtslibIndexError(
                "Format not indexable (htslib exit code -3)".to_string(),
            )),
            -4 => Err(DependencyError::HtslibIndexError(
                "Failed to create and/or save the index (htslib exit code -4)"
                    .to_string(),
            )),
            i => Err(DependencyError::HtslibIndexError(format!(
                "Unknown htslib exit code ({}) received",
                i
            ))),
        }
    }
}

/// https://github.com/samtools/htslib/blob/3c6f83f11e6e99cdcbed3541ae7d8d5f7786578d/htslib/faidx.h#L101
pub fn index_fasta(path: &Path) -> Result<PathBuf, DependencyError> {
    let fname = std::ffi::CString::new(path.to_path_buf().into_os_string().into_vec())
        .map_err(|_| {
            DependencyError::HtslibIndexError(
                "Failed to convert fasta path into a CString".to_string(),
            )
        })?;
    unsafe {
        match rust_htslib::htslib::fai_build(fname.as_ptr()) {
            0 => Ok(path.add_extension(OsStr::new(".fai"))),
            -1 => Err(DependencyError::HtslibIndexError(format!(
                "Fasta indexing failed (htslib exit code -1) {:?}",
                path
            ))),
            i => Err(DependencyError::HtslibIndexError(format!(
                "Unknown htslib exit code ({}) received",
                i
            ))),
        }
    }
}

pub fn dependency_dir() -> PathBuf {
    std::env::current_exe()
        .unwrap()
        .parent()
        .unwrap()
        .join("../../src/ext")
}

/// Checks whether the program is executable. If it is, it returns the full path to the
/// executable file
pub fn is_executable(program: &str) -> Option<String> {
    let cmd = format!("realpath $(command -v {})", program);
    let result = Command::new("sh").args(["-c", &cmd]).output();
    match result {
        Ok(output) => {
            let abspath = output.stdout.trim().to_str_lossy().to_string();
            if abspath.is_empty() {
                None
            } else {
                Some(abspath)
            }
        }
        _ => None,
    }
}

pub trait PathExt {
    fn add_extension(&self, extension: &OsStr) -> PathBuf;
    fn prefix(&self) -> Option<&str>;
}

impl PathExt for Path {
    fn add_extension(&self, extension: &OsStr) -> PathBuf {
        let mut s = self.as_os_str().to_os_string();
        s.push(extension);
        PathBuf::from(s)
    }
    /// Extracts the prefix (non-extension(s)) portion of [`self.file_name`]. This is a "left"
    /// variant of `file_stem` - meaning it takes the portion of the file name before the *first* `.`
    ///
    /// The prefix is:
    ///
    /// * [`None`], if there is no file name;
    /// * The entire file name if there is no embedded `.`;
    /// * The entire file name if the file name begins with `.` and has no other `.`s within;
    /// * Otherwise, the portion of the file name before the first `.`
    fn prefix(&self) -> Option<&str> {
        let stem = self.file_stem()?;
        let s = stem.to_str()?;
        if let Some(i) = s.find('.') {
            Some(&s[..i])
        } else {
            Some(s)
        }
    }
}

// Some extension methods for VCF Records as some of this functionality hasn't been released yet
pub trait VcfExt {
    fn end(&self) -> i64;
    fn rlen(&self) -> i64;
    fn range(&self) -> Range<i64>;
    fn contig(&self) -> String;
    fn coverage(&self) -> Option<(Vec<i32>, Vec<i32>)>;
    fn fraction_read_support(&self) -> Option<f32>;
    fn gt_conf(&self) -> Option<f32>;
    fn called_allele(&self) -> i32;
    fn is_pass(&self) -> bool;
    fn slice(&self, iv: &Range<i64>, ix: Option<usize>) -> &[u8];
    fn argmatch(&self, other: &Self) -> Option<usize>;
    fn is_indel(&self) -> bool;
}

impl VcfExt for bcf::Record {
    fn end(&self) -> i64 {
        self.pos() + self.rlen()
    }
    fn rlen(&self) -> i64 {
        self.inner().rlen
    }
    fn range(&self) -> Range<i64> {
        self.pos()..self.end()
    }

    fn contig(&self) -> String {
        String::from_utf8_lossy(&Vec::from(
            self.header()
                .rid2name(self.rid().expect("rid not set"))
                .expect("unable to find rid in header"),
        ))
        .to_string()
    }

    /// Returns the coverage on the forward and reverse strand (for the first sample only)
    fn coverage(&self) -> Option<(Vec<i32>, Vec<i32>)> {
        let fwd_covgs = self.format(Tags::FwdCovg.value()).integer().ok()?;
        let rev_covgs = self.format(Tags::RevCovg.value()).integer().ok()?;
        // we are making an assumption that we only ever deal with VCFs with one sample
        Some((fwd_covgs[0].to_owned(), rev_covgs[0].to_owned()))
    }

    fn fraction_read_support(&self) -> Option<f32> {
        let (fc, rc) = self.coverage()?;
        if fc.len() < 2 {
            return Some(1.0);
        }
        let gt = match self.called_allele() {
            i if i < 0 => return None,
            i => i,
        } as usize;
        let called_covg = (fc[gt] + rc[gt]) as f32;
        let mut other_covg = 0;

        if gt > 0 {
            other_covg = fc[0] + rc[0];
        } else {
            for (i, (f_cov, r_cov)) in fc.iter().zip(&rc).enumerate() {
                if i == gt {
                    continue;
                } else {
                    let cov: i32 = f_cov + r_cov;
                    if cov > other_covg {
                        other_covg = cov;
                    }
                }
            }
        }

        match called_covg / (called_covg + other_covg as f32) {
            f if f.is_nan() => None,
            f => Some(f),
        }
    }

    fn gt_conf(&self) -> Option<f32> {
        let gt_conf = self.format(Tags::GtypeConf.value()).float().ok()?;
        // there can only be one value for GT_CONF
        Some(gt_conf[0][0])
    }

    fn called_allele(&self) -> i32 {
        match self.genotypes() {
            Err(_) => -1,
            Ok(gts) => match gts.get(0)[..] {
                [GenotypeAllele::Unphased(i)] | [GenotypeAllele::Phased(i)] => i,
                _ => -1,
            },
        }
    }

    fn is_pass(&self) -> bool {
        self.has_filter(&Id(0))
    }

    /// Slice the specified allele with the given interval. If `None` is provided for the index,
    /// the called allele is used. If the called allele is NULL then REF is used.
    /// If `ix` is out of bounds for the number of alleles, an empty slice is returned.
    fn slice(&self, iv: &Range<i64>, ix: Option<usize>) -> &[u8] {
        let gt = match ix {
            None => match self.called_allele() {
                i if i < 0 => 0, // we assume REF for NULL calls
                i => i as usize,
            },
            Some(i) if i < self.allele_count() as usize => i,
            _ => return b"",
        };
        let allele = self.alleles()[gt];
        let allele_iv = self.pos()..self.pos() + allele.len() as i64;
        let isec = match allele_iv.intersect(iv) {
            Some(i) => {
                let s = (i.start - self.pos()) as usize;
                let e = min(s + (i.end - i.start) as usize, allele.len());
                s..e
            }
            None => return b"",
        };
        &self.alleles()[gt][isec]
    }

    /// Looks to match sequences between the record and another. This will match the called
    /// allele for Self with any of the alleles of `other`. Returns `Some` if there is a match, or
    /// `None` otherwise. If there are multiple matches, it will return the longest length match. If
    /// there are multiple matches of the same length, it will return the match with the biggest
    /// difference to the reference allele i.e., it will return the longest indel match
    fn argmatch(&self, other: &Self) -> Option<usize> {
        let called_len = match self.called_allele() {
            0 => self.rlen(),
            i if i > 0 => self.alleles()[i as usize].len() as i64,
            _ => return None,
        };
        let called_diff = (called_len - self.rlen()).abs();

        let mut match_ix = None;
        let mut match_diff: Option<i64> = None;

        let other_iv = self.pos()..(self.pos() + called_len);
        let other_ref = other.slice(&(self.pos()..i64::MAX), Some(0));
        let other_alleles = other.alleles();
        for (i, al) in other_alleles.iter().enumerate() {
            let is_indel = al.len() != other.alleles()[0].len();
            // we only want to compare snps with snp and indels with indels
            if self.is_indel() != is_indel {
                continue;
            }

            let iv = other.pos()..(other.pos() + al.len() as i64);
            let seq = self.slice(&iv, None);
            if seq.is_empty() {
                continue;
            }

            let other_seq = other.slice(&other_iv, Some(i));
            let diff = (other_ref.len() as i64 - al.len() as i64).abs() as i64;

            if seq != other_seq {
                continue;
            }

            if self.called_allele() == 0 && i == 0 {
                // the called allele is ref and we have a match with other's ref
                // short circuit and return a match with the ref - i.e. not resistant
                return Some(0);
            }
            if !self.is_indel() && !is_indel {
                // i.e. both S/MNPs
                let overlap_iv =
                    max(self.pos(), other.pos())..min(self.end(), other.end());
                let right_overhang_iv = overlap_iv.end..max(self.end(), other.end());
                let left_overhang_iv = min(self.pos(), other.pos())..overlap_iv.start;

                let self_overlap = self.slice(&overlap_iv, None);
                let self_left_overhang = if self.pos() == left_overhang_iv.start {
                    self.slice(&left_overhang_iv, None)
                } else {
                    other.slice(&left_overhang_iv, Some(0))
                };
                let self_right_overhang = if self.end() == right_overhang_iv.end {
                    self.slice(&right_overhang_iv, None)
                } else {
                    other.slice(&right_overhang_iv, Some(0))
                };
                let other_overlap = other.slice(&overlap_iv, Some(i));
                let other_left_overhang = if other.pos() == left_overhang_iv.start {
                    other.slice(&left_overhang_iv, Some(i))
                } else {
                    self.slice(&left_overhang_iv, Some(0))
                };
                let other_right_overhang = if other.end() == right_overhang_iv.end {
                    other.slice(&right_overhang_iv, Some(i))
                } else {
                    self.slice(&right_overhang_iv, Some(0))
                };

                let combined_self =
                    [self_left_overhang, self_overlap, self_right_overhang].concat();
                let combined_other =
                    [other_left_overhang, other_overlap, other_right_overhang].concat();

                if combined_other != combined_self {
                    continue;
                }
            }

            let diff_diff = (called_diff - diff).abs();
            match match_diff {
                Some(i) if i <= diff_diff => {}
                _ => {
                    match_diff = Some(diff_diff);
                    match_ix = Some(i);
                }
            }
        }
        match_ix
    }

    /// Tests if the record encodes an insertion/deletion (indel). This is only relevant if the
    /// record has called an alternate allele. Otherwise this function always returns `false`.
    fn is_indel(&self) -> bool {
        if self.called_allele() < 1 {
            false
        } else {
            self.alleles()[0].len()
                != self.alleles()[self.called_allele() as usize].len()
        }
    }
}

pub trait GffExt {
    fn name(&self) -> Option<&str>;
}

impl GffExt for gff::Record {
    fn name(&self) -> Option<&str> {
        match self
            .attributes()
            .iter()
            .find(|&entry| entry.key() == "Name")
        {
            None => None,
            Some(entry) => Some(entry.value()),
        }
    }
}

/// Reverse complement a DNA sequence
///
/// # Example
///
/// ```rust
/// use drprg::revcomp;
/// let seq = b"ATGCTTCCAGAA";
///
/// let actual = revcomp(seq);
/// let expected = b"TTCTGGAAGCAT";///
/// assert_eq!(actual, expected)
/// ```
///
/// # Note
/// Implementation is taken from https://doi.org/10.1101/082214
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|c| if c & 2 != 0 { c ^ 4 } else { c ^ 21 })
        .collect()
}

pub fn find_prg_index_in(dir: &Path) -> Option<PathBuf> {
    let mut fname: Option<String> = None;
    for entry in fs::read_dir(dir).ok()? {
        let path = entry.ok()?.path();
        if path.extension() == Some(OsStr::new("idx")) {
            fname = Some(path.file_name().unwrap().to_string_lossy().to_string());
        }
    }
    fname.map(|f| dir.join(f))
}

#[cfg(test)]
mod tests {
    use rust_htslib::bcf;
    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::Header;
    use std::any::Any;
    use std::io::Write;
    use tempfile::NamedTempFile;

    use crate::filter::test::{
        bcf_record_set_covg, bcf_record_set_gt, populate_bcf_header,
    };

    use super::*;

    #[test]
    fn test_revcomp() {
        let seq = b"ATGCTTCCAGAA";

        let actual = revcomp(seq);
        let expected = b"TTCTGGAAGCAT";

        assert_eq!(actual, expected)
    }

    #[test]
    fn path_is_executable() {
        let program = "ls";
        let executable = is_executable(program).unwrap();

        assert!(Path::new(&executable).is_absolute())
    }

    #[test]
    fn path_is_not_executable() {
        let program = "foobar";
        assert!(is_executable(program).is_none())
    }

    #[test]
    fn path_is_executable_resolves_full_path() {
        let program = "/bin/mkdir";
        let expected = Path::new(program)
            .canonicalize()
            .unwrap()
            .to_string_lossy()
            .to_string();
        assert_eq!(is_executable(program), Some(expected))
    }

    #[test]
    fn from_path_or_where_path_is_executable() {
        let path = Some(PathBuf::from("/bin/ls"));
        let default = PathBuf::from("DEFAULT");
        let actual = from_path_or(&path, &default).unwrap();

        let path = Path::new(&actual);
        assert!(path.is_absolute());
        assert_eq!(path.file_name().unwrap(), "ls")
    }

    #[test]
    fn from_path_or_where_path_isnt_executable_and_neither_is_default() {
        let path = Some(PathBuf::from("/bin/XZY"));
        let default = PathBuf::from("DEFAULT");
        let actual = from_path_or(&path, &default);
        assert!(actual.is_none())
    }

    #[test]
    fn from_path_or_where_path_is_none_and_default_isnt_executable() {
        let path = None;
        let default = PathBuf::from("DEFAULT");
        let actual = from_path_or(&path, &default);
        assert!(actual.is_none())
    }

    #[test]
    fn from_path_or_where_path_is_none_and_default_is_executable() {
        let path = None;
        let default = PathBuf::from("/bin/ls");
        let actual = from_path_or(&path, &default).unwrap();
        let path = Path::new(&actual);
        assert!(path.is_absolute());
        assert_eq!(path.file_name().unwrap(), "ls")
    }

    #[test]
    fn from_path_or_where_path_is_none_default_isnt_executable_but_default_filename_is()
    {
        let path = None;
        let default = PathBuf::from("/XZY/ls");
        let actual = from_path_or(&path, &default).unwrap();
        let path = Path::new(&actual);
        assert!(path.is_absolute());
        assert_eq!(path.file_name().unwrap(), "ls")
    }

    #[test]
    fn add_extension_empty_input() {
        let path = Path::new("foo.bar");

        let actual = path.add_extension("".as_ref());

        assert_eq!(actual, path)
    }

    #[test]
    fn add_extension_one_extension() {
        let path = Path::new("foo.bar");

        let actual = path.add_extension(".baz".as_ref());
        let expected = PathBuf::from("foo.bar.baz");

        assert_eq!(actual, expected)
    }

    #[test]
    fn add_extension_no_extension() {
        let path = Path::new("foo");

        let actual = path.add_extension(".baz".as_ref());
        let expected = PathBuf::from("foo.baz");

        assert_eq!(actual, expected)
    }

    #[test]
    fn file_prefix_no_file_ext() {
        let path = Path::new("dir/foo");

        let actual = path.prefix().unwrap();
        let expected = OsStr::new("foo");

        assert_eq!(actual, expected)
    }

    #[test]
    fn file_prefix_one_ext() {
        let path = Path::new("dir/foo.txt");

        let actual = path.prefix().unwrap();
        let expected = OsStr::new("foo");

        assert_eq!(actual, expected)
    }

    #[test]
    fn file_prefix_two_ext() {
        let path = Path::new("dir/foo.tar.gz");

        let actual = path.prefix().unwrap();
        let expected = OsStr::new("foo");

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_rlen() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        assert_eq!(record.rlen(), 0);
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        assert_eq!(record.rlen(), 3)
    }

    #[test]
    fn test_record_end() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record.set_pos(5);

        assert_eq!(record.end(), 8)
    }

    #[test]
    fn test_record_range() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record.set_pos(5);

        assert_eq!(record.range(), 5..8)
    }

    #[test]
    fn test_record_coverage() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[5, 0])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[6, 1])
            .expect("Failed to set reverse coverage");

        let actual = record.coverage();
        let expected = Some((vec![5, 0], vec![6, 1]));

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_coverage_no_tag() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample");
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");

        let actual = record.coverage();
        assert!(actual.is_none())
    }

    #[test]
    fn test_record_called_allele() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        let actual = record.called_allele();
        let expected = 1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_called_allele_is_null() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::UnphasedMissing])
            .unwrap();
        let actual = record.called_allele();
        let expected = -1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_gt_conf() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence">"#);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        record
            .push_format_float(Tags::GtypeConf.value(), &[3.4])
            .expect("Failed to set GT_CONF");

        let actual = record.gt_conf();
        let expected = Some(3.4);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_gt_conf_no_tag() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample");
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let record = vcf.empty_record();
        let actual = record.gt_conf();
        assert!(actual.is_none())
    }

    #[test]
    fn test_record_fraction_read_support() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[5, 0], &[4, 1]);
        bcf_record_set_gt(&mut record, 0);

        let actual = record.fraction_read_support();
        let expected = Some(0.9);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_fraction_read_support_alt() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[5, 0], &[4, 1]);
        bcf_record_set_gt(&mut record, 1);

        let actual = record.fraction_read_support();
        let expected = Some(0.1);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_fraction_read_support_zero_coverage() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[0, 0], &[0, 0]);
        bcf_record_set_gt(&mut record, 1);

        let actual = record.fraction_read_support();
        assert!(actual.is_none())
    }

    #[test]
    fn test_record_fraction_read_support_is_null() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[4, 4], &[0, 10]);
        bcf_record_set_gt(&mut record, -1);

        let actual = record.fraction_read_support();
        assert!(actual.is_none())
    }

    #[test]
    fn test_record_fraction_read_support_called_alt_compares_to_ref() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[4, 4, 7], &[0, 10, 1]);
        bcf_record_set_gt(&mut record, 1);

        let actual = record.fraction_read_support();
        let expected = Some(14.0 / (14.0 + 4.0));
        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_fraction_read_support_called_ref_compares_to_highest_alt() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[4, 4, 7], &[0, 10, 1]);
        bcf_record_set_gt(&mut record, 0);

        let actual = record.fraction_read_support();
        let expected = Some(4.0 / (14.0 + 4.0));
        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_is_pass() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        header
            .push_record(br#"##FILTER=<ID=foo,Description="sample is a foo fighter">"#);
        header.push_record(br#"##FILTER=<ID=bar,Description="a horse walks into...">"#);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        assert!(record.is_pass());
        record
            .push_filter(&record.header().name_to_id(b"foo").unwrap())
            .unwrap();
        assert!(!record.is_pass());
    }

    #[test]
    fn test_unwrap_or_continue() {
        let v = &[Ok(1), Ok(2), Err(3), Ok(4)];
        let mut sum = 0;
        let mut i = 0;
        for x in v {
            i += 1;
            sum += unwrap_or_continue!(x);
        }
        assert_eq!(i, 4);
        assert_eq!(sum, 7)
    }

    #[test]
    fn test_record_slice_ref_first_base() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        let iv = 0..1;

        let actual = record.slice(&iv, None);
        let expected = b"A";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_ref_last_base() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(0);
        let iv = 2..10;

        let actual = record.slice(&iv, None);
        let expected = b"G";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_ref_iv_spans_whole_and_more() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let iv = 2..10;

        let actual = record.slice(&iv, None);
        let expected = b"AGG";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_alt_iv_spans_whole_and_more() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(5);
        let iv = 2..10;

        let actual = record.slice(&iv, None);
        let expected = b"TG";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_alt_middle_base() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(5);
        let iv = 7..8;

        let actual = record.slice(&iv, None);
        let expected = b"A";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_empty_iv() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(5);
        let iv = 7..7;

        let actual = record.slice(&iv, None);
        let expected = b"";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_null_gt_uses_ref() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::UnphasedMissing])
            .unwrap();
        record.set_pos(5);
        let iv = 7..9;

        let actual = record.slice(&iv, None);
        let expected = b"G";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_no_iv_overlap_left() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let iv = 0..5;

        let actual = record.slice(&iv, None);
        let expected = b"";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_no_iv_overlap_right() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let iv = 8..10;

        let actual = record.slice(&iv, None);
        let expected = b"";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_specify_non_called_allele() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let iv = 6..110;
        let ix = Some(1);

        let actual = record.slice(&iv, ix);
        let expected = b"GAAA";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_specify_non_called_allele_out_of_bounds() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let iv = 6..110;
        let ix = Some(10);

        let actual = record.slice(&iv, ix);
        let expected = b"";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_specify_non_called_allele_mixed_lengths() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGAAG", b"TGAAAGGAAA", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let iv = 6..110;
        let ix = Some(2);

        let actual = record.slice(&iv, ix);
        let expected = b"";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_slice_specify_non_called_allele_single_base_olap() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ACG", b"AGAAA", b"GAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(7);
        let iv = 5..8;
        let ix = Some(2);

        let actual = record.slice(&iv, ix);
        let expected = b"G";

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_is_indel_no_alt() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();

        assert!(!record.is_indel())
    }

    #[test]
    fn test_record_is_indel_null() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"A"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(-1)])
            .unwrap();

        assert!(!record.is_indel())
    }

    #[test]
    fn test_record_is_indel_but_ref_call() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"A"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();

        assert!(!record.is_indel())
    }

    #[test]
    fn test_record_is_indel_deletion() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"A"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();

        assert!(record.is_indel())
    }

    #[test]
    fn test_record_is_indel_insertion() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"AAAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();

        assert!(record.is_indel())
    }

    #[test]
    fn test_record_is_indel_snp() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();

        assert!(!record.is_indel())
    }

    #[test]
    fn test_record_is_indel_mnp() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();

        assert!(!record.is_indel())
    }

    #[test]
    fn test_record_is_indel_snp_and_indel_in_alleles_called_snp() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"GA", b"G"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(2)])
            .unwrap();

        assert!(!record.is_indel())
    }

    #[test]
    fn test_record_is_indel_snp_and_indel_in_alleles_called_indel() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"GA", b"G"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();

        assert!(record.is_indel())
    }

    #[test]
    fn test_record_argmatch_same_record() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(5);

        let actual = record.argmatch(&other);
        let expected = Some(0);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_no_match() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ACG", b"AGAAA"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(5);

        let actual = record.argmatch(&other);
        let expected = None;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_only_overlap_matches() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ACG", b"AGAAA", b"G"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(7);

        let actual = record.argmatch(&other);
        let expected = None;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_self_spans_other() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"CCGGCATAT", b"CTGGCATAT", b"CCGGCAAAT"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(237);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"CGG", b"TGG"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(238);

        let actual = record.argmatch(&other);
        let expected = Some(1);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_only_overlap_matches_the_rest_doesnt() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ACG", b"AGAAA", b"GAAA"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(7);

        let actual = record.argmatch(&other);
        let expected = None;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_multiple_matches_at_overlap_returns_ref() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"GAA", b"GGAAA", b"GAAA"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(7);

        let actual = record.argmatch(&other);
        let expected = Some(0);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_multiple_matches_at_overlap_returns_shortest() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TTGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(5);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"GAAA", b"GGAAA", b"GA"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(7);

        let actual = record.argmatch(&other);
        let expected = Some(2);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_no_overlap() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TGAAA"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(5);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ACG", b"GGAAA", b"GAAA"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(9);

        let actual = record.argmatch(&other);
        let expected = None;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_single_base_deletion() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ATC", b"ACT", b"ACC", b"ACA", b"ACG", b"AC"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(5)])
            .unwrap();
        record.set_pos(161);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ATC", b"AC", b"ACC", b"ACA", b"ACG"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(161);

        let actual = record.argmatch(&other);
        let expected = Some(1);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_deletion_matches_longest() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"CCCCC", b"C"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(161);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"CCCCC", b"CCC", b"C"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(161);

        let actual = record.argmatch(&other);
        let expected = Some(2);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_deletion_matches_closest() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"CCCCC", b"CCC"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(161);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"CCCCC", b"CCCC", b"C"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(161);

        let actual = record.argmatch(&other);
        let expected = Some(1);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_deletion_matches_closest_overlap() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"CCCCC", b"CCC"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(160);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"CCCCC", b"CCCC", b"C"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(161);

        let actual = record.argmatch(&other);
        let expected = Some(1);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_single_base_insertion() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"G", b"GT"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(2197);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"G", b"GT", b"ACC", b"ACA", b"ACG"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(2197);

        let actual = record.argmatch(&other);
        let expected = Some(1);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_record_is_ref_and_matches_both_return_ref() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"ATTC"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record.set_pos(1396);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] =
            &[b"CTGAGCCAATTCATGGACCAGAACAACCC", b"CTGAGCCAACAGAACAACCC"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(1388);

        let actual = record.argmatch(&other);
        let expected = Some(0);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_insertion_matches_longest() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"C", b"CCCCC"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(161);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"C", b"CCC", b"CCCCC"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(161);

        let actual = record.argmatch(&other);
        let expected = Some(2);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_null_returns_none() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ATC", b"ACT", b"ACC", b"ACA", b"ACG", b"AC"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::UnphasedMissing])
            .unwrap();
        record.set_pos(161);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ATC", b"ACT", b"ACC", b"ACA", b"ACG"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(161);

        let actual = record.argmatch(&other);
        let expected = None;

        assert_eq!(actual, expected)
    }

    #[test]
    /// https://github.com/mbhall88/drprg/issues/11#issue-1312343270
    fn test_record_argmatch_overlap_base_matches_but_not_same() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ACGACG", b"ACGACA", b"GCGACG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(2)])
            .unwrap();
        record.set_pos(714);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"GCA", b"GAA", b"GAG"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(712);

        let actual = record.argmatch(&other);
        let expected = None;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_argmatch_overlap_base_matches_and_so_do_backfilled_seqs() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"ACGACG", b"ACGACA", b"GCGACG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(2)])
            .unwrap();
        record.set_pos(714);
        let mut other = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"GCA", b"GAA", b"GCG"];
        other.set_alleles(alleles).expect("Failed to set alleles");
        other
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        other.set_pos(712);

        let actual = record.argmatch(&other);
        let expected = Some(2);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_gff_record_name_is_present() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();

        assert_eq!(record.name(), Some("dnaA"))
    }

    #[test]
    fn test_gff_record_name_is_not_present() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF);
        let record = reader.records().next().unwrap().unwrap();

        assert!(record.name().is_none())
    }

    #[test]
    fn test_list_prgs_with_novel_variants() {
        let contents = r"1 samples
Sample ERR2510154
2 loci with denovo variants
gid
49 nodes
(0 [0, 89) AGTAAGCGATGCGTGGCCGAGCGGCTGGGCCAGCGTCTCGAGAGCGGAGAATGTTTCACGTGAAACATGACACAGACCTCACGAGCCGG)
(1 [92, 93) C)
(3 [100, 126) GGAGTGCGTAATGTCTCCGATCGAGC)
(4 [129, 130) C)
(6 [137, 143) CGCGGC)
(8 [174, 198) GTCTGCGATCTTCGGACCGCGGCG)
(10 [206, 218) TGGCCTTGCTCG)
(11 [222, 224) GC)
1 denovo variants for this locus
279     A
ahpC
9 nodes
(0 [0, 10) TAAATATGGT)
(1 [13, 16) GTG)
(4 [31, 36) ATATA)
(5 [39, 50) TCACCTTTGCC)
(9 [86, 106) TGACAGCGACTTCACGGCAC)
(10 [109, 110) G)
(12 [118, 393) ATGGAATGTCGCAACCAAATGCATTGTCCGCTTTGATGATGAGGAGAGTCATGCCACTGCTAACCATTGGCGATCAATTCCCCGCCTACCAGCTCACCGCTCTCATCGGCGGTGACCTGTCCAAGGTCGACGCCAAGCAGCCCGGCGACTACTTCACCACTATCACCAGTGACGAACACCCAGGCAAGTGGCGGGTGGTGTTCTTTTGGCCGAAAGACTTCACGTTCGTGTGCCCTACCGAGATCGCGGCGTTCAGCAAGCTCAATGACGAGTTC)
(13 [397, 398) G)
(15 [407, 869) AGGACCGCGACGCCCAGATCCTGGGGGTTTCGATTGACAGCGAATTCGCGCATTTCCAGTGGCGTGCACAGCACAACGACCTCAAAACGTTACCCTTCCCGATGCTCTCCGACATCAAGCGCGAACTCAGCCAAGCCGCAGGTGTCCTCAACGCCGACGGTGTGGCCGACCGCGTGACCTTTATCGTCGACCCCAACAACGAGATCCAGTTCGTCTCGGCCACCGCCGGTTCGGTGGGACGCAACGTCGATGAGGTACTGCGAGTGCTCGACGCCCTCCAGTCCGACGAGCTGTGCGCATGCAACTGGCGCAAGGGCGACCCGACGCTAGACGCTGGCGAACTCCTCAAGGCTTCGGCCTAACCGGGATCTGGTTGGCCGGGAATCAATGAGTATAGAAAAGCTCAAGGCCGCGCTCCCCGAGTACGCCAAAGACATCAAGCTGAACCTGAGCTCAATCACC)
2 denovo variants for this locus
246     T       G
249     T       G".to_string();
        let mut tmpfile = NamedTempFile::new().unwrap();
        let path = tmpfile.path().to_owned();
        {
            tmpfile.write_all(contents.as_ref()).unwrap();
        }

        let actual = Pandora::list_prgs_with_novel_variants(&path).unwrap();
        let expected = vec!["gid".to_string(), "ahpC".to_string()];

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_list_prgs_with_novel_variants_different_num_vars() {
        let contents = r"1 samples
Sample ERR2510154
4 loci with denovo variants
gid
49 nodes
(0 [0, 89) AGTAAGCGATGCGTGGCCGAGCGGCTGGGCCAGCGTCTCGAGAGCGGAGAATGTTTCACGTGAAACATGACACAGACCTCACGAGCCGG)
(1 [92, 93) C)
(3 [100, 126) GGAGTGCGTAATGTCTCCGATCGAGC)
(4 [129, 130) C)
(6 [137, 143) CGCGGC)
(8 [174, 198) GTCTGCGATCTTCGGACCGCGGCG)
(10 [206, 218) TGGCCTTGCTCG)
(11 [222, 224) GC)
1 denovo variants for this locus
279     A
ahpC
9 nodes
(0 [0, 10) TAAATATGGT)
(1 [13, 16) GTG)
(4 [31, 36) ATATA)
(5 [39, 50) TCACCTTTGCC)
(9 [86, 106) TGACAGCGACTTCACGGCAC)
(10 [109, 110) G)
(12 [118, 393) ATGGAATGTCGCAACCAAATGCATTGTCCGCTTTGATGATGAGGAGAGTCATGCCACTGCTAACCATTGGCGATCAATTCCCCGCCTACCAGCTCACCGCTCTCATCGGCGGTGACCTGTCCAAGGTCGACGCCAAGCAGCCCGGCGACTACTTCACCACTATCACCAGTGACGAACACCCAGGCAAGTGGCGGGTGGTGTTCTTTTGGCCGAAAGACTTCACGTTCGTGTGCCCTACCGAGATCGCGGCGTTCAGCAAGCTCAATGACGAGTTC)
(13 [397, 398) G)
(15 [407, 869) AGGACCGCGACGCCCAGATCCTGGGGGTTTCGATTGACAGCGAATTCGCGCATTTCCAGTGGCGTGCACAGCACAACGACCTCAAAACGTTACCCTTCCCGATGCTCTCCGACATCAAGCGCGAACTCAGCCAAGCCGCAGGTGTCCTCAACGCCGACGGTGTGGCCGACCGCGTGACCTTTATCGTCGACCCCAACAACGAGATCCAGTTCGTCTCGGCCACCGCCGGTTCGGTGGGACGCAACGTCGATGAGGTACTGCGAGTGCTCGACGCCCTCCAGTCCGACGAGCTGTGCGCATGCAACTGGCGCAAGGGCGACCCGACGCTAGACGCTGGCGAACTCCTCAAGGCTTCGGCCTAACCGGGATCTGGTTGGCCGGGAATCAATGAGTATAGAAAAGCTCAAGGCCGCGCTCCCCGAGTACGCCAAAGACATCAAGCTGAACCTGAGCTCAATCACC)
2 denovo variants for this locus
246     T       G
249     T       G".to_string();
        let mut tmpfile = NamedTempFile::new().unwrap();
        let path = tmpfile.path().to_owned();
        {
            tmpfile.write_all(contents.as_ref()).unwrap();
        }

        let actual = Pandora::list_prgs_with_novel_variants(&path).unwrap_err();
        let expected = DependencyError::NovelVariantParsingError("".to_string());

        assert_eq!(actual.type_id(), expected.type_id())
    }
}
