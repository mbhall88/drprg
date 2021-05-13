use crate::VcfExt;
use float_cmp::approx_eq;
use rust_htslib::bcf;
use rust_htslib::bcf::Record;
use std::str::FromStr;
use thiserror::Error;

/// A collection of custom errors relating to the Tags enum
#[derive(Error, Debug, PartialEq)]
pub enum TagsError {
    /// The string is not known
    #[error("Unknown tag given {0}")]
    UnknownTag(String),
}

/// A collection of known VCF tags/fields/filters
#[derive(Debug, Ord, PartialOrd, Eq, PartialEq)]
pub enum Tags {
    FwdCovg,
    RevCovg,
    LowCovg,
    HighCovg,
    StrandBias,
    Gaps,
    HighGaps,
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
            Tags::FwdCovg => b"MED_FWD_COVG",
            Tags::RevCovg => b"MED_REV_COVG",
            Tags::LowCovg => b"ld",
            Tags::HighCovg => b"hd",
            Tags::StrandBias => b"sb",
            Tags::Gaps => b"GAPS",
            Tags::HighGaps => b"hg",
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
    type Err = TagsError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "MED_FWD_COVG" => Ok(Tags::FwdCovg),
            "MED_REV_COVG" => Ok(Tags::RevCovg),
            "ld" => Ok(Tags::LowCovg),
            "hd" => Ok(Tags::HighCovg),
            "sb" => Ok(Tags::StrandBias),
            "GAPS" => Ok(Tags::Gaps),
            "hg" => Ok(Tags::HighGaps),
            "GT_CONF" => Ok(Tags::GtypeConf),
            "lgc" => Ok(Tags::LowGtConf),
            "lindel" => Ok(Tags::LongIndel),
            "frs" => Ok(Tags::LowSupport),
            "PASS" => Ok(Tags::Pass),
            "FT" => Ok(Tags::FormatFilter),
            "FAIL" => Ok(Tags::AllFail),
            _ => Err(TagsError::UnknownTag(s.to_string())),
        }
    }
}

#[derive(Default)]
struct FilterStatus {
    low_covg: bool,
    high_covg: bool,
    strand_bias: bool,
    low_gt_conf: bool,
    high_gaps: bool,
    long_indel: bool,
    low_support: bool,
}

impl FilterStatus {
    fn new() -> Self {
        FilterStatus::default()
    }
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
        if self.high_gaps {
            tags.push(Tags::HighGaps);
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
}

/// A trait to allow for filtering of a VCF. The provided method is `filter`
pub trait Filter {
    fn is_low_covg(&self, record: &bcf::Record) -> bool;
    fn is_high_covg(&self, record: &bcf::Record) -> bool;
    fn is_low_gt_conf(&self, record: &bcf::Record) -> bool;
    fn is_low_support(&self, record: &bcf::Record) -> bool;
    fn is_long_indel(&self, record: &bcf::Record) -> bool;
    fn filter(&self, record: &mut bcf::Record) {
        todo!()
    }
}

/// A struct to hold all of the values for the filter
pub struct Filterer {
    min_covg: i32,
    max_covg: i32,
    min_strand_bias: f32,
    min_gt_conf: f32,
    max_indel: Option<i32>,
    min_support: f32,
}

impl Default for Filterer {
    fn default() -> Self {
        Filterer {
            min_covg: -1,
            max_covg: i32::MAX,
            min_strand_bias: -1.0,
            min_gt_conf: -1.0,
            max_indel: None,
            min_support: -1.0,
        }
    }
}

impl Filter for Filterer {
    fn is_low_covg(&self, record: &Record) -> bool {
        // safe to unwrap as we are providing a default
        let (fc, rc) = record.coverage().unwrap_or_else(|| (vec![0], vec![0]));
        let default_gt = 0; // i.e. if GT is null, we use ref coverage
        let gt = match record.called_allele() {
            i if i < 0 => default_gt,
            i => i,
        } as usize;
        let covg = fc.get(gt).unwrap_or(&0) + rc.get(gt).unwrap_or(&0);
        covg < self.min_covg
    }

    fn is_high_covg(&self, record: &Record) -> bool {
        // safe to unwrap as we are providing a default
        let (fc, rc) = record.coverage().unwrap_or_else(|| (vec![0], vec![0]));
        let default_gt = 0; // i.e. if GT is null, we use ref coverage
        let gt = match record.called_allele() {
            i if i < 0 => default_gt,
            i => i,
        } as usize;
        let covg = fc.get(gt).unwrap_or(&0) + rc.get(gt).unwrap_or(&0);
        covg > self.max_covg
    }

    fn is_low_gt_conf(&self, record: &Record) -> bool {
        // we assume that a missing GT_CONF value implies 0
        let gt_conf = record.gt_conf().unwrap_or(0.0);
        gt_conf < self.min_gt_conf && !approx_eq!(f32, gt_conf, self.min_gt_conf)
    }

    fn is_low_support(&self, record: &Record) -> bool {
        match record.fraction_read_support() {
            Some(frs) => {
                frs < self.min_support && !approx_eq!(f32, frs, self.min_support)
            }
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
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;
    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::Header;
    use tempfile::NamedTempFile;

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
        let expected = Err(TagsError::UnknownTag(s.to_string()));

        assert_eq!(actual, expected)
    }

    #[test]
    fn filter_status_tags_none_fail() {
        let status = FilterStatus::new();

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
            high_gaps: true,
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
            Tags::HighGaps,
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
            high_gaps: true,
            long_indel: false,
            low_support: false,
        };

        let mut actual = status.tags();
        actual.sort_unstable();
        let mut expected = vec![Tags::StrandBias, Tags::HighGaps];
        expected.sort_unstable();

        assert_eq!(actual, expected)
    }

    pub(crate) fn populate_bcf_header(header: &mut bcf::Header) {
        header.push_sample(b"sample")
            .push_record(br#"##FORMAT=<ID=MED_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#)
            .push_record(br#"##FORMAT=<ID=MED_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        ).push_record(br#"##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence">"#);
    }

    pub(crate) fn bcf_record_set_covg(
        record: &mut bcf::Record,
        fwd_covg: &[i32],
        rev_covg: &[i32],
    ) {
        let alleles: &[&[u8]] = &[b"AGG", b"TG"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_format_integer(b"MED_FWD_COVG", fwd_covg)
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MED_REV_COVG", rev_covg)
            .expect("Failed to set reverse coverage");
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

        let mut filt = Filterer::default();
        filt.min_covg = 2;

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

        let mut filt = Filterer::default();
        filt.min_covg = 2;

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

        let mut filt = Filterer::default();
        filt.min_covg = 3;

        assert!(filt.is_low_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_low_covg_null_gt_ref_below() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 3], &[1, 3]);
        bcf_record_set_gt(&mut record, -1);

        let mut filt = Filterer::default();
        filt.min_covg = 3;

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

        let mut filt = Filterer::default();
        filt.min_covg = 3;

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

        let mut filt = Filterer::default();
        filt.min_covg = 3;

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

        let mut filt = Filterer::default();
        filt.max_covg = 2;

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

        let mut filt = Filterer::default();
        filt.max_covg = 2;

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

        let mut filt = Filterer::default();
        filt.max_covg = 3;

        assert!(filt.is_high_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_high_covg_null_gt_ref_above() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[10, 3], &[1, 3]);
        bcf_record_set_gt(&mut record, -1);

        let mut filt = Filterer::default();
        filt.max_covg = 3;

        assert!(filt.is_high_covg(&record))
    }

    #[test]
    fn filter_bcf_record_is_high_covg_null_gt_ref_below() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();
        populate_bcf_header(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, &[1, 3], &[1, 3]);
        bcf_record_set_gt(&mut record, -1);

        let mut filt = Filterer::default();
        filt.max_covg = 3;

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

        let mut filt = Filterer::default();
        filt.max_covg = 3;

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

        let mut filt = Filterer::default();
        filt.min_gt_conf = 3.2;

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

        let mut filt = Filterer::default();
        filt.min_gt_conf = 31.2;

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

        let mut filt = Filterer::default();
        filt.min_gt_conf = 11.1;

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

        let mut filt = Filterer::default();
        filt.min_support = 0.9;

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

        let mut filt = Filterer::default();
        filt.min_support = 0.9;

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

        let mut filt = Filterer::default();
        filt.min_support = 0.90;

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

        let mut filt = Filterer::default();
        filt.min_support = 0.90;

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

        let mut filt = Filterer::default();
        filt.min_support = 0.90;

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

        let mut filt = Filterer::default();
        filt.max_indel = Some(1);

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

        let mut filt = Filterer::default();
        filt.max_indel = Some(1);

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

        let mut filt = Filterer::default();
        filt.max_indel = Some(10);

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

        let mut filt = Filterer::default();
        filt.max_indel = Some(10);

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

        let mut filt = Filterer::default();
        filt.max_indel = Some(2);

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

        let mut filt = Filterer::default();
        filt.max_indel = Some(2);

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

        let mut filt = Filterer::default();
        filt.max_indel = Some(0);

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

        let mut filt = Filterer::default();
        filt.max_indel = Some(0);

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

        let mut filt = Filterer::default();
        filt.max_indel = Some(0);

        assert!(!filt.is_long_indel(&record))
    }
}
