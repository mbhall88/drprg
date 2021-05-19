use std::ffi::OsStr;
use std::fs::File;
use std::ops::Range;
use std::os::unix::ffi::OsStringExt;
use std::path::{Path, PathBuf};
use std::process::Command;

use bstr::ByteSlice;
use log::{debug, error};
use rust_htslib::bcf;
use rust_htslib::bcf::header::Id;
use rust_htslib::bcf::record::GenotypeAllele;
use thiserror::Error;

use crate::filter::Tags;

pub mod filter;
mod interval;

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
    /// An error associated with indexing VCFs with htslib
    #[error("Failed to index VCF with htslib: {0}")]
    HtslibIndexError(String),
}

pub struct Bcftools {
    executable: String,
}

impl Bcftools {
    pub fn from_path(path: &Option<PathBuf>) -> Result<Self, DependencyError> {
        let default = dependency_dir().join(BCFTOOLS_BIN);
        let executable = from_path_or(&path, &default);
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
}

pub struct MakePrg {
    executable: String,
}

impl MakePrg {
    /// Creates a MakePrg object from a path or from default if None given
    pub fn from_path(path: &Option<PathBuf>) -> Result<MakePrg, DependencyError> {
        let default = dependency_dir().join(MAKE_PRG_BIN);
        let executable = from_path_or(&path, &default);
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
        output_update_ds: &Path,
        args: I,
    ) -> Result<(), DependencyError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let dir = tempfile::tempdir().map_err(DependencyError::ProcessError)?;
        let prefix = "dr";
        let cmd_result = Command::new(&self.executable)
            .current_dir(&dir)
            .arg("from_msa")
            .args(args)
            .args(&["-o", prefix, "-i"])
            .arg(input)
            .output();

        match cmd_result {
            Ok(cmd_output) if !cmd_output.status.success() => {
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
                let tmp_update_ds = tmp_prefix.with_extension("update_DS");

                // have to use fs::copy here as fs::rename fails inside a container as the tempdir is
                // not on the same "mount" as the local filesystem see more info at
                // https://doc.rust-lang.org/std/fs/fn.rename.html#platform-specific-behavior
                std::fs::copy(tmp_prg, output_prg)
                    .map_err(|source| DependencyError::FileError { source })?;
                std::fs::copy(tmp_update_ds, output_update_ds)
                    .map_err(|source| DependencyError::FileError { source })?;

                Ok(())
            }
            Err(err) => {
                error!("make_prg failed to run with error: {}", err.to_string());
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
        let executable = from_path_or(&path, &default);
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
        let cmd_output = Command::new(&self.executable)
            .arg("index")
            .args(args)
            .arg(input)
            .output()
            .map_err(DependencyError::ProcessError)?;

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
        let fixed_args = &[
            "map",
            "--genotype",
            "-o",
            &outdir.to_string_lossy(),
            "-g",
            &MTB_GENOME_SIZE.to_string(),
            "--vcf-refs",
            &vcf_ref.to_string_lossy(),
        ];
        let cmd_output = Command::new(&self.executable)
            .args(fixed_args)
            .args(args)
            .arg(prg)
            .arg(reads)
            .output()
            .map_err(DependencyError::ProcessError)?;

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
}

pub struct MultipleSeqAligner {
    executable: String,
}

impl MultipleSeqAligner {
    pub fn from_path(
        path: &Option<PathBuf>,
    ) -> Result<MultipleSeqAligner, DependencyError> {
        let default = dependency_dir().join(MAFFT_BIN);
        let executable = from_path_or(&path, &default);
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
        let result = Command::new(&self.executable)
            .args(args)
            .arg(input)
            .stdout(ostream)
            .output();
        match result {
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
                .or_else(|| is_executable(&*default_fname))
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
    let result = Command::new("sh").args(&["-c", &cmd]).output();
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
    fn file_prefix(&self) -> Option<&str>;
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
    fn file_prefix(&self) -> Option<&str> {
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
    fn coverage(&self) -> Option<(Vec<i32>, Vec<i32>)>;
    fn fraction_read_support(&self) -> Option<f32>;
    fn gt_conf(&self) -> Option<f32>;
    fn called_allele(&self) -> i32;
    fn is_pass(&self) -> bool;
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

    /// Returns the coverage on the forward and reverse strand (for the first sample only)
    fn coverage(&self) -> Option<(Vec<i32>, Vec<i32>)> {
        let fwd_covgs = self.format(Tags::FwdCovg.value()).integer().ok()?;
        let rev_covgs = self.format(Tags::RevCovg.value()).integer().ok()?;
        // we are making an assumption that we only ever deal with VCFs with one sample
        Some((fwd_covgs[0].to_owned(), rev_covgs[0].to_owned()))
    }

    fn fraction_read_support(&self) -> Option<f32> {
        let (fc, rc) = self.coverage()?;
        let total_covg = fc.iter().chain(&rc).sum::<i32>() as f32;
        let gt = match self.called_allele() {
            i if i < 0 => return None,
            i => i,
        } as usize;
        let called_covg = (fc[gt] + rc[gt]) as f32;
        match called_covg / total_covg {
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
        self.has_filter(Id(0))
    }
}

#[cfg(test)]
mod tests {
    use rust_htslib::bcf;
    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::Header;
    use tempfile::NamedTempFile;

    use crate::filter::test::{
        bcf_record_set_covg, bcf_record_set_gt, populate_bcf_header,
    };

    use super::*;

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
        let program = "src/ext/pandora";
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

        let actual = path.file_prefix().unwrap();
        let expected = OsStr::new("foo");

        assert_eq!(actual, expected)
    }

    #[test]
    fn file_prefix_one_ext() {
        let path = Path::new("dir/foo.txt");

        let actual = path.file_prefix().unwrap();
        let expected = OsStr::new("foo");

        assert_eq!(actual, expected)
    }

    #[test]
    fn file_prefix_two_ext() {
        let path = Path::new("dir/foo.tar.gz");

        let actual = path.file_prefix().unwrap();
        let expected = OsStr::new("foo");

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_record_rlen() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MED_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MED_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_format_integer(b"MED_FWD_COVG", &[5, 0])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MED_REV_COVG", &[6, 1])
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[4, 4], &[0, 10]);
        bcf_record_set_gt(&mut record, -1);

        let actual = record.fraction_read_support();
        assert!(actual.is_none())
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
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        assert!(record.is_pass());
        record.push_filter(record.header().name_to_id(b"foo").unwrap());
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
}
