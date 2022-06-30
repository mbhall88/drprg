use std::io::Write;
use std::str::FromStr;

use float_cmp::approx_eq;
use lazy_static::lazy_static;
use rust_htslib::bcf;
use rust_htslib::bcf::header::Id;
use rust_htslib::bcf::Record;
use structopt::StructOpt;
use thiserror::Error;

use crate::VcfExt;

const MIN_COVG: i32 = -1;
const MAX_COVG: i32 = i32::MAX;
const MIN_SB: f32 = -1.0;
const MIN_GTCONF: f32 = -1.0;
const MAX_INDEL: Option<i32> = None;
const MIN_FRS: f32 = -1.0;
lazy_static! {
    static ref MIN_COVG_STR: String = MIN_COVG.to_string();
    static ref MAX_COVG_STR: String = MAX_COVG.to_string();
    static ref MIN_SB_STR: String = MIN_SB.to_string();
    static ref MIN_GTCONF_STR: String = MIN_GTCONF.to_string();
    static ref MIN_FRS_STR: String = MIN_FRS.to_string();
}

/// A collection of custom errors relating to the Tags enum
#[derive(Error, Debug, PartialEq)]
pub enum FilterError {
    /// The string is not known
    #[error("Unknown tag given {0}")]
    UnknownTag(String),
    /// Tried to add a filter tag not in the header
    #[error("Tried to add filter {0}, which is not in the header")]
    TagNotInHeader(String),
}

/// A collection of known VCF tags/fields/filters
#[derive(Debug, Ord, PartialOrd, Eq, PartialEq)]
pub enum Tags {
    FwdCovg,
    RevCovg,
    LowCovg,
    HighCovg,
    StrandBias,
    GtypeConf,
    LowGtConf,
    LongIndel,
    LowSupport,
    Pass,
    FormatFilter,
    AllFail,
}

impl Tags {
    /// Returns the string representation of the tag
    pub fn value(&self) -> &[u8] {
        match *self {
            Tags::FwdCovg => b"MEAN_FWD_COVG",
            Tags::RevCovg => b"MEAN_REV_COVG",
            Tags::LowCovg => b"ld",
            Tags::HighCovg => b"hd",
            Tags::StrandBias => b"sb",
            Tags::GtypeConf => b"GT_CONF",
            Tags::LowGtConf => b"lgc",
            Tags::LongIndel => b"lindel",
            Tags::LowSupport => b"frs",
            Tags::Pass => b"PASS",
            Tags::FormatFilter => b"FT",
            Tags::AllFail => b"FAIL",
        }
    }
}

impl FromStr for Tags {
    type Err = FilterError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "MEAN_FWD_COVG" => Ok(Tags::FwdCovg),
            "MEAN_REV_COVG" => Ok(Tags::RevCovg),
            "ld" => Ok(Tags::LowCovg),
            "hd" => Ok(Tags::HighCovg),
            "sb" => Ok(Tags::StrandBias),
            "GT_CONF" => Ok(Tags::GtypeConf),
            "lgc" => Ok(Tags::LowGtConf),
            "lindel" => Ok(Tags::LongIndel),
            "frs" => Ok(Tags::LowSupport),
            "PASS" => Ok(Tags::Pass),
            "FT" => Ok(Tags::FormatFilter),
            "FAIL" => Ok(Tags::AllFail),
            _ => Err(FilterError::UnknownTag(s.to_string())),
        }
    }
}

#[derive(Default)]
struct FilterStatus {
    low_covg: bool,
    high_covg: bool,
    strand_bias: bool,
    low_gt_conf: bool,
    long_indel: bool,
    low_support: bool,
}

impl FilterStatus {
    /// Returns the Tags for the filters this object fails
    fn tags(&self) -> Vec<Tags> {
        let mut tags = vec![];
        if self.low_covg {
            tags.push(Tags::LowCovg);
        }
        if self.high_covg {
            tags.push(Tags::HighCovg);
        }
        if self.low_gt_conf {
            tags.push(Tags::LowGtConf);
        }
        if self.strand_bias {
            tags.push(Tags::StrandBias);
        }
        if self.long_indel {
            tags.push(Tags::LongIndel);
        }
        if self.low_support {
            tags.push(Tags::LowSupport);
        }
        if tags.is_empty() {
            tags.push(Tags::Pass);
        }
        tags
    }

    fn apply_filters_to(&self, record: &mut bcf::Record) -> Result<(), FilterError> {
        let mut ids: Vec<Id> = vec![];
        for tag in self.tags() {
            let id = record.header().name_to_id(tag.value()).map_err(|_| {
                FilterError::TagNotInHeader(
                    String::from_utf8_lossy(tag.value()).to_string(),
                )
            })?;
            ids.push(id);
        }
        record.set_filters(&ids);
        Ok(())
    }
}

/// A trait to allow for filtering of a VCF. The provided method is `filter`
pub trait Filter {
    fn _covg_for_gt(&self, record: &Record) -> i32;
    fn is_low_covg(&self, record: &bcf::Record) -> bool;
    fn is_high_covg(&self, record: &bcf::Record) -> bool;
    fn is_low_gt_conf(&self, record: &bcf::Record) -> bool;
    fn is_low_support(&self, record: &bcf::Record) -> bool;
    fn is_long_indel(&self, record: &bcf::Record) -> bool;
    fn has_strand_bias(&self, record: &bcf::Record) -> bool;
    fn filter(&self, record: &mut bcf::Record) -> Result<(), FilterError> {
        let status = FilterStatus {
            low_covg: self.is_low_covg(record),
            high_covg: self.is_high_covg(record),
            strand_bias: self.has_strand_bias(record),
            low_gt_conf: self.is_low_gt_conf(record),
            long_indel: self.is_long_indel(record),
            low_support: self.is_low_support(record),
        };

        status.apply_filters_to(record)?;
        Ok(())
    }
}

/// A struct to hold all of the values for the filter
#[derive(StructOpt, Debug)]
pub struct Filterer {
    /// Minimum depth of coverage allowed on variants
    #[structopt(short = "d", long, default_value = &MIN_COVG_STR, hidden_short_help = true, value_name = "INT", allow_hyphen_values = true)]
    pub min_covg: i32,
    /// Maximum depth of coverage allowed on variants
    #[structopt(short = "D", long, default_value = &MAX_COVG_STR, hidden_short_help = true, value_name = "INT", allow_hyphen_values = true)]
    pub max_covg: i32,
    /// Minimum strand bias ratio allowed on variants
    ///
    /// For example, setting to 0.25 requires >=25% of total (allele) coverage on both
    /// strands for an allele.
    #[structopt(short = "b", long, default_value = &MIN_SB_STR, hidden_short_help = true, value_name = "FLOAT", allow_hyphen_values = true)]
    pub min_strand_bias: f32,
    /// Minimum genotype confidence (GT_CONF) score allow on variants
    #[structopt(short = "g", long, default_value = &MIN_GTCONF_STR, hidden_short_help = true, value_name = "FLOAT, allow_hyphen_values = true")]
    pub min_gt_conf: f32,
    /// Maximum (absolute) length of insertions/deletions allowed
    #[structopt(
        short = "L",
        long,
        hidden_short_help = true,
        value_name = "INT",
        allow_hyphen_values = true
    )]
    pub max_indel: Option<i32>,
    /// Minimum fraction of read support
    ///
    /// For example, setting to 0.9 requires >=90% of coverage for the variant to be on the called
    /// allele
    #[structopt(short = "K", long, default_value = &MIN_FRS_STR, hidden_short_help = true, value_name = "FLOAT", allow_hyphen_values = true)]
    pub min_frs: f32,
}

impl Default for Filterer {
    fn default() -> Self {
        Filterer {
            min_covg: MIN_COVG,
            max_covg: MAX_COVG,
            min_strand_bias: MIN_SB,
            min_gt_conf: MIN_GTCONF,
            max_indel: MAX_INDEL,
            min_frs: MIN_FRS,
        }
    }
}

impl Filter for Filterer {
    fn _covg_for_gt(&self, record: &Record) -> i32 {
        // safe to unwrap as we are providing a default
        let (fc, rc) = record.coverage().unwrap_or_else(|| (vec![0], vec![0]));
        let gt = record.called_allele();
        let covg: i32 = if gt < 0 {
            // if GT is null, we use total covg on all alleles
            fc.iter().chain(rc.iter()).sum()
        } else {
            fc.get(gt as usize).unwrap_or(&0) + rc.get(gt as usize).unwrap_or(&0)
        };
        covg
    }

    fn is_low_covg(&self, record: &Record) -> bool {
        let covg = self._covg_for_gt(record);
        covg < self.min_covg
    }

    fn is_high_covg(&self, record: &Record) -> bool {
        let covg = self._covg_for_gt(record);
        covg > self.max_covg
    }

    fn is_low_gt_conf(&self, record: &Record) -> bool {
        // we assume that a missing GT_CONF value implies 0
        let gt_conf = record.gt_conf().unwrap_or(0.0);
        gt_conf < self.min_gt_conf && !approx_eq!(f32, gt_conf, self.min_gt_conf)
    }

    fn is_low_support(&self, record: &Record) -> bool {
        match record.fraction_read_support() {
            Some(frs) => frs < self.min_frs && !approx_eq!(f32, frs, self.min_frs),
            None => false,
        }
    }

    fn is_long_indel(&self, record: &Record) -> bool {
        let gt = record.called_allele();
        if gt < 1 || self.max_indel.is_none() {
            false
        } else {
            let l = match record.alleles().get(gt as usize) {
                Some(a) => a.len(),
                None => 0,
            } as i64;

            let indel_len = (record.rlen() - l).abs();
            indel_len > self.max_indel.unwrap() as i64
        }
    }

    fn has_strand_bias(&self, record: &Record) -> bool {
        type B = Vec<f32>;
        let (fc, rc): (B, B) = match record.coverage() {
            Some((f, r)) => (
                f.iter().map(|i| *i as f32).collect::<B>(),
                r.iter().map(|i| *i as f32).collect::<B>(),
            ),
            None => return false,
        };
        let ratio = match record.called_allele() {
            -1 => {
                let total_fc: f32 = fc.iter().sum();
                let total_rc: f32 = rc.iter().sum();
                let total = total_fc + total_rc;

                if approx_eq!(f32, total, 0.0) {
                    None
                } else {
                    Some(total_fc.min(total_rc) / total)
                }
            }
            gt => {
                let sum_covg = fc[gt as usize] + rc[gt as usize];
                if approx_eq!(f32, sum_covg, 0.0) {
                    None
                } else {
                    Some(fc[gt as usize].min(rc[gt as usize]) / sum_covg)
                }
            }
        };
        match ratio {
            Some(r) => {
                r < self.min_strand_bias && !approx_eq!(f32, r, self.min_strand_bias)
            }
            None => false,
        }
    }
}

impl Filterer {
    fn meta_info_line(id: &[u8], description: &[u8]) -> Vec<u8> {
        let mut line: Vec<u8> = b"##FILTER=<ID=".to_vec();
        line.write_all(id).unwrap();
        line.write_all(br#",Description=""#).unwrap();
        line.write_all(description).unwrap();
        line.write_all(br#"">"#).unwrap();
        line
    }
    /// Write meta-information lines to a VCF for the filters that are currently enabled
    pub fn add_filter_headers(&self, header: &mut bcf::Header) {
        if self.min_covg > MIN_COVG {
            let id = Tags::LowCovg.value();
            let desc =
                format!("Kmer coverage on called allele less than {}", self.min_covg);
            header.push_record(Self::meta_info_line(id, desc.as_bytes()).as_slice());
        }
        if self.max_covg < MAX_COVG {
            let id = Tags::HighCovg.value();
            let desc =
                format!("Kmer coverage on called allele more than {}", self.min_covg);
            header.push_record(Self::meta_info_line(id, desc.as_bytes()).as_slice());
        }
        if self.min_strand_bias > MIN_SB {
            let id = Tags::StrandBias.value();
            let desc = format!(
                "A strand on the called allele has less than {:.2}% of the coverage for that allele", self.min_strand_bias * 100.0,
            );
            header.push_record(Self::meta_info_line(id, desc.as_bytes()).as_slice());
        }
        if self.min_gt_conf > MIN_GTCONF {
            let id = Tags::LowGtConf.value();
            let desc = format!(
                "Genotype confidence score less than {:.1}",
                self.min_gt_conf,
            );
            header.push_record(Self::meta_info_line(id, desc.as_bytes()).as_slice());
        }
        if self.max_indel.is_some() {
            let id = Tags::LongIndel.value();
            let desc = format!("Indel is longer than {}bp", self.max_indel.unwrap(),);
            header.push_record(Self::meta_info_line(id, desc.as_bytes()).as_slice());
        }
        if self.min_frs > MIN_FRS {
            let id = Tags::LowSupport.value();
            let desc = format!(
                "Fraction of read support on called allele is less than {}",
                self.min_frs
            );
            header.push_record(Self::meta_info_line(id, desc.as_bytes()).as_slice());
        }
    }
}

#[cfg(test)]
pub(crate) mod test {
    use bstr::ByteSlice;
    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::Header;
    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn tags_value() {
        assert_eq!(Tags::GtypeConf.value(), b"GT_CONF")
    }

    #[test]
    fn tags_from_str() {
        let s = "frs";

        let actual = Tags::from_str(s);
        let expected = Ok(Tags::LowSupport);

        assert_eq!(actual, expected)
    }

    #[test]
    fn tags_from_str_unknown_tag() {
        let s = "foo";

        let actual = Tags::from_str(s);
        let expected = Err(FilterError::UnknownTag(s.to_string()));

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_status_tags_none_fail() {
        let status = FilterStatus::default();

        let actual = status.tags();
        let expected = &[Tags::Pass];

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_status_tags_all_fail() {
        let status = FilterStatus {
            low_covg: true,
            high_covg: true,
            strand_bias: true,
            low_gt_conf: true,
            long_indel: true,
            low_support: true,
        };

        let mut actual = status.tags();
        actual.sort_unstable();
        let mut expected = vec![
            Tags::LowCovg,
            Tags::HighCovg,
            Tags::LowGtConf,
            Tags::StrandBias,
            Tags::LongIndel,
            Tags::LowSupport,
        ];
        expected.sort_unstable();

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_status_tags_some_fail() {
        let status = FilterStatus {
            low_covg: false,
            high_covg: false,
            strand_bias: true,
            low_gt_conf: false,
            long_indel: false,
            low_support: true,
        };

        let mut actual = status.tags();
        actual.sort_unstable();
        let mut expected = vec![Tags::StrandBias, Tags::LowSupport];
        expected.sort_unstable();

        assert_eq!(actual, expected)
    }

    pub(crate) fn populate_bcf_header(header: &mut bcf::Header) {
        header.push_sample(b"sample")
            .push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#)
            .push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        ).push_record(br#"##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence">"#);
        header.push_record(br#"##contig=<ID=chrom>"#);
    }

    pub(crate) fn bcf_record_set_covg(
        record: &mut bcf::Record,
        fwd_covg: &[i32],
        rev_covg: &[i32],
    ) {
        record
            .push_format_integer(b"MEAN_FWD_COVG", fwd_covg)
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", rev_covg)
            .expect("Failed to set reverse coverage");
    }

    pub(crate) fn bcf_record_set_alleles(record: &mut bcf::Record, alleles: &[&[u8]]) {
        record.set_alleles(alleles).expect("Failed to set alleles");
    }

    pub(crate) fn bcf_record_set_gt(record: &mut bcf::Record, gt: i32) {
        let gts = match gt {
            i if i < 0 => vec![GenotypeAllele::UnphasedMissing],
            i => vec![GenotypeAllele::Unphased(i)],
        };
        record.push_genotypes(&gts).unwrap();
    }

    pub(crate) fn bcf_record_set_gt_conf(record: &mut bcf::Record, gt_conf: f32) {
        record
            .push_format_float(Tags::GtypeConf.value(), &[gt_conf])
            .expect("Failed to set gt conf");
    }

    #[test]
    fn filter_bcf_record_is_low_covg() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[5], &[5]);
        bcf_record_set_gt(&mut record, 0);

        let filt = Filterer {
            min_covg: 2,
            ..Default::default()
        };

        assert!(!filt.is_low_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_covg_same_as_min() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1], &[1]);
        bcf_record_set_gt(&mut record, 0);

        let filt = Filterer {
            min_covg: 2,
            ..Default::default()
        };

        assert!(!filt.is_low_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_covg_below_min() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1], &[1]);
        bcf_record_set_gt(&mut record, 0);

        let filt = Filterer {
            min_covg: 3,
            ..Default::default()
        };

        assert!(filt.is_low_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_covg_null_gt_sum_of_all_is_below() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 3], &[1, 3]);
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer {
            min_covg: 9,
            ..Default::default()
        };

        assert!(filt.is_low_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_covg_null_gt_ref_above() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[6, 3], &[1, 3]);
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer {
            min_covg: 3,
            ..Default::default()
        };

        assert!(!filt.is_low_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_covg_no_covg() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        header.remove_format(Tags::FwdCovg.value());
        header.remove_format(Tags::RevCovg.value());
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer {
            min_covg: 3,
            ..Default::default()
        };

        assert!(filt.is_low_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_covg_no_covg_but_min_covg_unset() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        header.remove_format(Tags::FwdCovg.value());
        header.remove_format(Tags::RevCovg.value());
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer::default();

        assert!(!filt.is_low_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_high_covg() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[5], &[5]);
        bcf_record_set_gt(&mut record, 0);

        let filt = Filterer {
            max_covg: 2,
            ..Default::default()
        };

        assert!(filt.is_high_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_high_covg_same_as_max() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1], &[1]);
        bcf_record_set_gt(&mut record, 0);

        let filt = Filterer {
            max_covg: 2,
            ..Default::default()
        };

        assert!(!filt.is_high_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_high_covg_above_max() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[2], &[2]);
        bcf_record_set_gt(&mut record, 0);

        let filt = Filterer {
            max_covg: 3,
            ..Default::default()
        };

        assert!(filt.is_high_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_high_covg_null_gt_sum_is_above() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[10, 3], &[1, 3]);
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer {
            max_covg: 12,
            ..Default::default()
        };

        assert!(filt.is_high_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_high_covg_null_gt_sum_of_all_is_above() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 3], &[1, 3]);
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer {
            max_covg: 9,
            ..Default::default()
        };

        assert!(!filt.is_high_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_high_covg_no_covg() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        header.remove_format(Tags::FwdCovg.value());
        header.remove_format(Tags::RevCovg.value());
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer {
            max_covg: 3,
            ..Default::default()
        };

        assert!(!filt.is_high_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_high_covg_no_covg_but_max_covg_unset() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        header.remove_format(Tags::FwdCovg.value());
        header.remove_format(Tags::RevCovg.value());
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer::default();

        assert!(!filt.is_high_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_gt_conf_no_gt_conf_but_gt_conf_unset() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        header.remove_format(Tags::GtypeConf.value());
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let record = vcf.empty_record();

        let filt = Filterer::default();

        assert!(!filt.is_low_gt_conf(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_gt_conf_above() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let conf = 11.1;
        bcf_record_set_gt_conf(&mut record, conf);

        let filt = Filterer {
            min_gt_conf: 3.2,
            ..Default::default()
        };

        assert!(!filt.is_low_gt_conf(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_gt_conf_below() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let conf = 11.1;
        bcf_record_set_gt_conf(&mut record, conf);

        let filt = Filterer {
            min_gt_conf: 31.2,
            ..Default::default()
        };

        assert!(filt.is_low_gt_conf(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_gt_conf_same_as_min() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let conf = 11.1;
        bcf_record_set_gt_conf(&mut record, conf);

        let filt = Filterer {
            min_gt_conf: 11.1,
            ..Default::default()
        };

        assert!(!filt.is_low_gt_conf(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_gt_conf_unset() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let conf = 0.0;
        bcf_record_set_gt_conf(&mut record, conf);

        let filt = Filterer::default();

        assert!(!filt.is_low_gt_conf(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_support_above() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 30], &[1, 31]);
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer {
            min_frs: 0.90,
            ..Default::default()
        };

        assert!(!filt.is_low_support(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_support_below() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 30], &[1, 31]);
        bcf_record_set_gt(&mut record, 0);

        let filt = Filterer {
            min_frs: 0.90,
            ..Default::default()
        };

        assert!(filt.is_low_support(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_support_same_as_min() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[5, 0], &[4, 1]);
        bcf_record_set_gt(&mut record, 0);

        let filt = Filterer {
            min_frs: 0.90,
            ..Default::default()
        };

        assert!(!filt.is_low_support(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_support_null_gt() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[5, 10], &[4, 1]);
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer {
            min_frs: 0.90,
            ..Default::default()
        };

        assert!(!filt.is_low_support(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_support_no_coverage() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[0, 0], &[0, 0]);
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer {
            min_frs: 0.90,
            ..Default::default()
        };

        assert!(!filt.is_low_support(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_support_no_min_set() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[5, 5], &[5, 5]);
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer::default();

        assert!(!filt.is_low_support(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_ins_above() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TCGAG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer {
            max_indel: Some(1),
            ..Default::default()
        };

        assert!(filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_del_above() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer {
            max_indel: Some(1),
            ..Default::default()
        };

        assert!(filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_ins_below() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TCGAG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer {
            max_indel: Some(10),
            ..Default::default()
        };

        assert!(!filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_del_below() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer {
            max_indel: Some(10),
            ..Default::default()
        };

        assert!(!filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_ins_same_as_max() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"TCGAG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer {
            max_indel: Some(2),
            ..Default::default()
        };

        assert!(!filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_del_below_same_as_max() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer {
            max_indel: Some(2),
            ..Default::default()
        };

        assert!(!filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_ref() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 0);

        let filt = Filterer {
            max_indel: Some(0),
            ..Default::default()
        };

        assert!(!filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_null() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, -1);

        let filt = Filterer {
            max_indel: Some(0),
            ..Default::default()
        };

        assert!(!filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_unset() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"AGG", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer::default();

        assert!(!filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_is_long_indel_snp() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);

        let filt = Filterer {
            max_indel: Some(0),
            ..Default::default()
        };

        assert!(!filt.is_long_indel(&record))
    }

    #[test]
    fn filter_bcf_record_has_strand_bias_above() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 30], &[1, 31]);
        bcf_record_set_gt(&mut record, 0);
        bcf_record_set_alleles(&mut record, &[b"A", b"T"]);

        let filt = Filterer {
            min_strand_bias: 0.1,
            ..Default::default()
        };

        assert!(!filt.has_strand_bias(&record))
    }

    #[test]
    fn filter_bcf_record_has_strand_bias_below() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 1], &[1, 9]);
        bcf_record_set_gt(&mut record, 1);
        bcf_record_set_alleles(&mut record, &[b"A", b"T"]);

        let filt = Filterer {
            min_strand_bias: 0.2,
            ..Default::default()
        };

        assert!(filt.has_strand_bias(&record))
    }

    #[test]
    fn filter_bcf_record_has_strand_bias_same_as_min() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 1], &[1, 9]);
        bcf_record_set_gt(&mut record, 1);
        bcf_record_set_alleles(&mut record, &[b"A", b"T"]);

        let filt = Filterer {
            min_strand_bias: 0.1,
            ..Default::default()
        };

        assert!(!filt.has_strand_bias(&record))
    }

    #[test]
    fn filter_bcf_record_has_strand_bias_no_coverage() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 0], &[1, 0]);
        bcf_record_set_gt(&mut record, 1);
        bcf_record_set_alleles(&mut record, &[b"A", b"T"]);

        let filt = Filterer {
            min_strand_bias: 0.1,
            ..Default::default()
        };

        assert!(!filt.has_strand_bias(&record))
    }

    #[test]
    fn filter_bcf_record_has_strand_bias_null_uses_all() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 0, 4], &[1, 0, 40]);
        bcf_record_set_gt(&mut record, -1);
        bcf_record_set_alleles(&mut record, &[b"A", b"T", b"C"]);

        let filt = Filterer {
            min_strand_bias: 0.1,
            ..Default::default()
        };

        assert!(!filt.has_strand_bias(&record));

        let filt2 = Filterer {
            min_strand_bias: 0.13,
            ..Default::default()
        };

        assert!(filt2.has_strand_bias(&record))
    }

    #[test]
    fn filter_all_pass() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        bcf_record_set_covg(&mut record, &[1, 0], &[1, 0]);
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);
        record.set_pos(1);

        let filt = Filterer::default();
        filt.filter(&mut record).unwrap();

        assert!(record.is_pass())
    }

    #[test]
    fn filter_low_covg() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        header.push_record(br#"##FILTER=<ID=ld,Description="low covg">"#);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        bcf_record_set_covg(&mut record, &[1, 10], &[1, 0]);
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 1);
        record.set_pos(1);

        let filt = Filterer {
            min_covg: 200,
            ..Default::default()
        };
        filt.filter(&mut record).unwrap();

        let id = record.header().name_to_id(Tags::LowCovg.value()).unwrap();
        assert!(record.has_filter(id));
        assert!(!record.is_pass())
    }

    #[test]
    fn filter_high_covg_and_low_support() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        header.push_record(br#"##FILTER=<ID=hd,Description="low covg">"#);
        header.push_record(br#"##FILTER=<ID=frs,Description="low covg">"#);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        bcf_record_set_covg(&mut record, &[10, 10], &[10, 0]);
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 0);
        record.set_pos(1);

        let filt = Filterer {
            max_covg: 2,
            min_frs: 0.95,
            ..Default::default()
        };
        filt.filter(&mut record).unwrap();

        let mut id = record.header().name_to_id(Tags::HighCovg.value()).unwrap();
        assert!(record.has_filter(id));
        id = record
            .header()
            .name_to_id(Tags::LowSupport.value())
            .unwrap();
        assert!(record.has_filter(id));
        assert!(!record.is_pass())
    }

    #[test]
    fn filter_status_apply_filters_to_missing_tag() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        header.push_record(br#"##FILTER=<ID=hd,Description="low covg">"#);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        bcf_record_set_covg(&mut record, &[10, 10], &[10, 0]);
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 0);
        record.set_pos(1);

        let stat = FilterStatus {
            low_covg: false,
            high_covg: true,
            strand_bias: false,
            low_gt_conf: false,
            long_indel: true,
            low_support: false,
        };

        let actual = stat.apply_filters_to(&mut record).unwrap_err();
        let expected = FilterError::TagNotInHeader(
            Tags::LongIndel.value().to_str_lossy().to_string(),
        );

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_status_apply_filters_to_tags_all_present() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        header.push_record(br#"##FILTER=<ID=sb,Description="low covg">"#);
        header.push_record(br#"##FILTER=<ID=lgc,Description="low covg">"#);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        bcf_record_set_covg(&mut record, &[10, 10], &[10, 0]);
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 0);
        record.set_pos(1);

        let stat = FilterStatus {
            low_covg: false,
            high_covg: false,
            strand_bias: true,
            low_gt_conf: true,
            long_indel: false,
            low_support: false,
        };

        stat.apply_filters_to(&mut record).unwrap();

        let mut id = record.header().name_to_id(Tags::LowGtConf.value()).unwrap();
        assert!(record.has_filter(id));
        id = record
            .header()
            .name_to_id(Tags::StrandBias.value())
            .unwrap();
        assert!(record.has_filter(id));
    }

    #[test]
    fn filter_status_apply_filters_to_all_pass() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        bcf_record_set_covg(&mut record, &[10, 10], &[10, 0]);
        record.set_alleles(alleles).expect("Failed to set alleles");
        bcf_record_set_gt(&mut record, 0);
        record.set_pos(1);

        let stat = FilterStatus::default();

        stat.apply_filters_to(&mut record).unwrap();

        assert!(record.is_pass())
    }

    #[test]
    fn filterer_add_filters_to_header_all_default_set_nothing() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        let filt = Filterer::default();
        filt.add_filter_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let view = vcf.header();

        assert!(view.name_to_id(Tags::LowGtConf.value()).is_err())
    }

    #[test]
    fn filterer_add_filters_to_header_all_set() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        let filt = Filterer {
            min_covg: 0,
            max_covg: 0,
            min_strand_bias: 0.0,
            min_gt_conf: 0.0,
            max_indel: Some(1),
            min_frs: 0.0,
        };
        filt.add_filter_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let view = vcf.header();

        assert!(view.name_to_id(Tags::LowGtConf.value()).is_ok());
        assert!(view.name_to_id(Tags::LowCovg.value()).is_ok());
        assert!(view.name_to_id(Tags::HighCovg.value()).is_ok());
        assert!(view.name_to_id(Tags::LowSupport.value()).is_ok());
        assert!(view.name_to_id(Tags::LongIndel.value()).is_ok());
        assert!(view.name_to_id(Tags::StrandBias.value()).is_ok())
    }

    #[test]
    fn filterer_add_filters_to_header_some_set() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        let filt = Filterer {
            min_covg: MIN_COVG,
            max_covg: 0,
            min_strand_bias: MIN_SB,
            min_gt_conf: 0.0,
            max_indel: Some(1),
            min_frs: 0.0,
        };
        filt.add_filter_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let view = vcf.header();

        assert!(view.name_to_id(Tags::LowGtConf.value()).is_ok());
        assert!(view.name_to_id(Tags::LowCovg.value()).is_err());
        assert!(view.name_to_id(Tags::HighCovg.value()).is_ok());
        assert!(view.name_to_id(Tags::LowSupport.value()).is_ok());
        assert!(view.name_to_id(Tags::LongIndel.value()).is_ok());
        assert!(view.name_to_id(Tags::StrandBias.value()).is_err())
    }
}
