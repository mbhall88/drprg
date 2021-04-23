use std::collections::{HashMap, HashSet};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::str::FromStr;

use regex::Regex;
use serde::{de, Deserialize, Deserializer};
use thiserror::Error;

pub(crate) type Panel = HashMap<String, HashSet<PanelRecord>>;

/// A collection of custom errors relating to the command line interface for this package.
#[derive(Error, Debug, PartialEq)]
pub enum PanelError {
    /// Unrecognised residue type
    #[error("{0} is not a known residue type")]
    UnknownResidueType(String),
    /// The variant string is not in the correct format
    #[error("The variant is not in the correct format [<STR><INT><STR>]: {0}")]
    IncorrectVariantFormat(String),
}

/// An enum representing the panel residue types we recognise
#[derive(Debug, PartialEq, Hash, Eq, Clone)]
pub enum Residue {
    Nucleic,
    Amino,
}

impl fmt::Display for Residue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match &self {
            Residue::Nucleic => "DNA",
            Residue::Amino => "PROT",
        };
        write!(f, "{}", s)
    }
}

impl FromStr for Residue {
    type Err = PanelError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "DNA" => Ok(Residue::Nucleic),
            "PROT" => Ok(Residue::Amino),
            unknown => Err(PanelError::UnknownResidueType(unknown.to_string())),
        }
    }
}

impl<'de> Deserialize<'de> for Residue {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        FromStr::from_str(&s).map_err(de::Error::custom)
    }
}

/// A variant object holding information about the reference, position, and new allele
#[derive(Debug, PartialEq, Hash, Eq)]
pub struct Variant {
    pub reference: String,
    pub pos: i64,
    pub new: String,
}

impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", &self.reference, self.pos, &self.new)
    }
}

impl FromStr for Variant {
    type Err = PanelError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        lazy_static! {
            static ref VARIANT_REGEX: Regex =
                Regex::new(r"^([a-zA-Z]+)(-?\d+)([a-zA-Z]+)$").unwrap();
        }
        let caps = VARIANT_REGEX.captures(s);
        match caps {
            // we use len==4 since every regex has at least one capture group - the full match
            Some(c) if c.len() == 4 => Ok(Variant {
                reference: c[1].to_string(),
                new: c[3].to_string(),
                pos: c.get(2).unwrap().as_str().parse::<i64>().unwrap(),
            }),
            _ => Err(PanelError::IncorrectVariantFormat(s.to_string())),
        }
    }
}

impl<'de> Deserialize<'de> for Variant {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        FromStr::from_str(&s).map_err(de::Error::custom)
    }
}

/// An object for working with panel mutations
#[derive(Debug)]
pub struct Mutation {
    residue: Residue,
    name: String,
    reference: String,
    new: String,
    pos: i64,
}

impl Mutation {
    /// Create a mutation from a [`PanelRecord`].
    pub fn from_record(record: &PanelRecord) -> Self {
        let name = format!("{}_{}", record.gene, record.variant);
        let residue = record.residue.to_owned();
        let reference = record.variant.reference.to_owned();
        let pos = record.variant.pos;
        let new = record.variant.new.to_owned();
        Mutation {
            residue,
            name,
            reference,
            new,
            pos,
        }
    }
}

/// An object containing the data in a row of the panel file
#[derive(Debug, Deserialize)]
pub struct PanelRecord {
    /// The gene the mutation is within
    pub gene: String,
    /// The variant describing reference, position, alternate (i.e. A4G)
    pub variant: Variant,
    /// The molecular alphabet the mutation uses - i.e. Nucleic or Amino Acid
    pub residue: Residue,
    /// The drugs the mutation causes resistance to
    #[serde(deserialize_with = "str_to_set")]
    pub drugs: HashSet<String>,
}

impl Hash for PanelRecord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.gene.hash(state);
        self.variant.hash(state);
        self.residue.hash(state);
    }
}

impl PartialEq for PanelRecord {
    fn eq(&self, other: &Self) -> bool {
        self.gene == other.gene
            && self.variant == other.variant
            && self.residue == other.residue
    }
}

impl Eq for PanelRecord {}

/// Allow serde to deserialize the drugs section of the panel
fn str_to_set<'de, D>(deserializer: D) -> Result<HashSet<String>, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    Ok(s.split(',').map(|item| item.to_owned()).collect())
}

#[cfg(test)]
mod tests {
    use std::iter::FromIterator;

    use super::*;

    #[test]
    fn mutation_type_display() {
        assert_eq!(Residue::Amino.to_string(), "PROT");
        assert_eq!(Residue::Nucleic.to_string(), "DNA")
    }

    #[test]
    fn mutation_type_from_str() {
        assert_eq!(Residue::from_str("DNA"), Ok(Residue::Nucleic));
        assert_eq!(Residue::from_str("PROT"), Ok(Residue::Amino));
        assert_eq!(Residue::from_str("prot"), Ok(Residue::Amino));
        assert!(Residue::from_str("foo").is_err())
    }

    #[test]
    fn variant_display() {
        let var = Variant {
            reference: "K".to_string(),
            pos: -6,
            new: "A".to_string(),
        };
        let actual = var.to_string();
        let expected = String::from("K-6A");
        assert_eq!(actual, expected)
    }

    #[test]
    fn variant_from_str_with_negative_pos() {
        let s = "K-1Q";
        let actual = Variant::from_str(s).unwrap();
        let expected = Variant {
            reference: "K".to_string(),
            pos: -1,
            new: "Q".to_string(),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn variant_from_str_with_positive_pos() {
        let s = "K11Q";
        let actual = Variant::from_str(s).unwrap();
        let expected = Variant {
            reference: "K".to_string(),
            pos: 11,
            new: "Q".to_string(),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn variant_from_str_with_multi_ref() {
        let s = "AT11C";
        let actual = Variant::from_str(s).unwrap();
        let expected = Variant {
            reference: "AT".to_string(),
            pos: 11,
            new: "C".to_string(),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn variant_from_str_with_multi_new() {
        let s = "AT11CGG";
        let actual = Variant::from_str(s).unwrap();
        let expected = Variant {
            reference: "AT".to_string(),
            pos: 11,
            new: "CGG".to_string(),
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn variant_from_str_with_missing_pos() {
        let s = "ATCGG";
        let result = Variant::from_str(s);
        assert!(result.is_err())
    }

    #[test]
    fn variant_from_str_with_missing_ref() {
        let s = "5ATCGG";
        let result = Variant::from_str(s);
        assert!(result.is_err())
    }

    #[test]
    fn variant_from_str_with_missing_new() {
        let s = "ATCGG1";
        let result = Variant::from_str(s);
        assert!(result.is_err())
    }

    #[test]
    fn variant_from_str_with_invalid_ref() {
        let s = "+6T";
        let result = Variant::from_str(s);
        assert!(result.is_err())
    }

    #[test]
    fn variant_from_str_with_invalid_new() {
        let s = "T6 ";
        let result = Variant::from_str(s);
        assert!(result.is_err())
    }

    #[test]
    fn deserialise_correct_panel_record() {
        const CSV: &[u8] = b"gene\tK1S\tPROT\td1,d2";

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(CSV);

        let record: PanelRecord = reader.deserialize().next().unwrap().unwrap();

        let expected = PanelRecord {
            gene: "gene".to_string(),
            variant: Variant::from_str("K1S").unwrap(),
            residue: Residue::Amino,
            drugs: HashSet::from_iter(vec!["d1".to_string(), "d2".to_string()]),
        };

        assert_eq!(record, expected)
    }

    #[test]
    fn deserialise_panel_record_with_unknown_mutation_type() {
        const CSV: &[u8] = b"gene\tK1S\tfoo\td1,d2";

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(CSV);

        let result: Result<PanelRecord, csv::Error> =
            reader.deserialize().next().unwrap();
        assert!(result.is_err())
    }

    #[test]
    fn deserialise_panel_record_with_wrong_delim() {
        const CSV: &[u8] = b"gene;K1S;foo;d1,d2";

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(CSV);

        let result: Result<PanelRecord, csv::Error> =
            reader.deserialize().next().unwrap();
        assert!(result.is_err())
    }

    #[test]
    fn deserialise_panel_record_with_missing_field() {
        const CSV: &[u8] = b"gene\tDNA\td1,d2";

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(CSV);

        let result: Result<PanelRecord, csv::Error> =
            reader.deserialize().next().unwrap();
        assert!(result.is_err())
    }
}
