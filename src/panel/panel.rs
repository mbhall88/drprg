use serde::{de, Deserialize, Deserializer};
use std::collections::HashSet;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::str::FromStr;
use thiserror::Error;

pub(crate) type Panel = HashSet<PanelRecord>;

/// A collection of custom errors relating to the command line interface for this package.
#[derive(Error, Debug, PartialEq)]
pub enum PanelError {
    /// Unrecognised mutation type
    #[error("{0} is not known mutation type")]
    UnknownMutationType(String),
}

/// An enum representing the panel mutation types we recognise
#[derive(Debug, PartialEq, Hash, Eq)]
pub enum MutationType {
    Nucleic,
    Amino,
}

impl fmt::Display for MutationType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match &self {
            MutationType::Nucleic => "DNA",
            MutationType::Amino => "PROT",
        };
        write!(f, "{}", s)
    }
}

impl FromStr for MutationType {
    type Err = PanelError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "DNA" => Ok(MutationType::Nucleic),
            "PROT" => Ok(MutationType::Amino),
            unknown => Err(PanelError::UnknownMutationType(unknown.to_string())),
        }
    }
}

impl<'de> Deserialize<'de> for MutationType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        FromStr::from_str(&s).map_err(de::Error::custom)
    }
}

/// An object containing the data in a row of the panel file
#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct PanelRecord {
    /// The gene the mutation is within
    pub gene: String,
    /// The mutation describing reference, position, alternate (i.e. A4G)
    pub mutation: String,
    /// The molecular alphabet the mutation uses - i.e. Nucleic or Amino Acid
    pub mutation_type: MutationType,
    /// The drugs the mutation causes resistance to
    #[serde(deserialize_with = "str_to_set")]
    pub drugs: HashSet<String>,
}

impl Hash for PanelRecord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.gene.hash(state);
        self.mutation.hash(state);
        self.mutation_type.hash(state);
    }
}

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
    use super::*;
    use std::iter::FromIterator;

    #[test]
    fn mutation_type_display() {
        assert_eq!(MutationType::Amino.to_string(), "PROT");
        assert_eq!(MutationType::Nucleic.to_string(), "DNA")
    }

    #[test]
    fn mutation_type_from_str() {
        assert_eq!(MutationType::from_str("DNA"), Ok(MutationType::Nucleic));
        assert_eq!(MutationType::from_str("PROT"), Ok(MutationType::Amino));
        assert_eq!(MutationType::from_str("prot"), Ok(MutationType::Amino));
        assert!(MutationType::from_str("foo").is_err())
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
            mutation: "K1S".to_string(),
            mutation_type: MutationType::Amino,
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
