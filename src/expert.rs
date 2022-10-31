use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
use std::collections::HashSet;
use std::fmt;
use std::str::FromStr;
use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum ExpertError {
    /// Unrecognised variant type
    #[error("{0} is not a recognised variant type")]
    UnknownVariantType(String),
}

enum VariantType {
    Frameshift,
    Nonsense,
    Missense,
}

impl fmt::Display for VariantType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match &self {
            VariantType::Frameshift => "frameshift",
            VariantType::Missense => "missense",
            VariantType::Nonsense => "nonsense",
        };
        write!(f, "{}", s)
    }
}

impl Serialize for VariantType {
    fn serialize<S>(
        &self,
        serializer: S,
    ) -> Result<<S as Serializer>::Ok, <S as Serializer>::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

impl FromStr for VariantType {
    type Err = ExpertError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "frameshift" => Ok(VariantType::Frameshift),
            "nonsense" => Ok(VariantType::Nonsense),
            "missense" => Ok(VariantType::Missense),
            unknown => Err(ExpertError::UnknownVariantType(unknown.to_string())),
        }
    }
}

impl<'de> Deserialize<'de> for VariantType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        FromStr::from_str(&s).map_err(de::Error::custom)
    }
}

pub struct Rule {
    variant_type: VariantType,
    gene: String,
    start: isize,
    end: isize,
    drugs: HashSet<String>,
}
