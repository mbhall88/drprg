use std::str::FromStr;
use thiserror::Error;

/// A collection of custom errors relating to the Tags enum
#[derive(Error, Debug, PartialEq)]
pub enum TagsError {
    /// The string is not known
    #[error("Unknown tag given {0}")]
    UnknownTag(String),
}
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

pub trait Filter {}

#[cfg(test)]
mod test {
    use super::*;

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
}
