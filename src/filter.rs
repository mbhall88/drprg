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
enum Tags {
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
    fn value(&self) -> &str {
        match *self {
            Tags::FwdCovg => "MED_FWD_COVG",
            Tags::RevCovg => "MED_REV_COVG",
            Tags::LowCovg => "ld",
            Tags::HighCovg => "hd",
            Tags::StrandBias => "sb",
            Tags::Gaps => "GAPS",
            Tags::HighGaps => "hg",
            Tags::GtypeConf => "GT_CONF",
            Tags::LowGtConf => "lgc",
            Tags::LongIndel => "lindel",
            Tags::LowSupport => "frs",
            Tags::Pass => "PASS",
            Tags::FormatFilter => "FT",
            Tags::AllFail => "FAIL",
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
        FilterStatus {
            low_covg: false,
            high_covg: false,
            strand_bias: false,
            low_gt_conf: false,
            high_gaps: false,
            long_indel: false,
            low_support: false,
        }
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
    fn is_high_gaps(&self, record: &bcf::Record) -> bool;
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

impl Filter for Filterer {
    fn is_low_covg(&self, record: &Record) -> bool {
        todo!()
    }

    fn is_high_covg(&self, record: &Record) -> bool {
        todo!()
    }

    fn is_low_gt_conf(&self, record: &Record) -> bool {
        todo!()
    }

    fn is_low_support(&self, record: &Record) -> bool {
        todo!()
    }

    fn is_high_gaps(&self, record: &Record) -> bool {
        todo!()
    }

    fn is_long_indel(&self, record: &Record) -> bool {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rust_htslib::bcf::Header;
    use tempfile::NamedTempFile;

    #[test]
    fn tags_value() {
        assert_eq!(Tags::GtypeConf.value(), "GT_CONF")
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

    fn bcf_record_set_covg(record: &mut bcf::Record, fwd_covg: i32, rev_covg: i32) {
        todo!()
    }

    #[test]
    fn filter_bcf_record_is_low_covg() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let header = Header::new();
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::VCF).unwrap();
        let mut record = vcf.empty_record();
        bcf_record_set_covg(&mut record, 5, 5);
    }
}
