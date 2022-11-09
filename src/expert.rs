use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
use serde_derive::Deserialize;
use std::collections::BTreeSet;
use std::fmt;
use std::hash::Hash;
use std::str::FromStr;
use thiserror::Error;

#[derive(Error, Debug, PartialEq, Eq)]
pub enum ExpertError {
    /// Unrecognised variant type
    #[error("{0} is not a recognised variant type")]
    UnknownVariantType(String),
}

#[derive(Eq, PartialEq, Debug, Hash)]
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

/// An object containing the data in a row of the expert rules file
#[derive(Debug, Deserialize, Eq, PartialEq, Hash)]
pub struct Rule {
    /// The variant type the rule describes
    variant_type: VariantType,
    /// The gene the rule applies to
    gene: String,
    /// An optional codon start coordinate in the gene to apply rule from (1-based inclusive).
    /// If not given, the start of the gene is inferred
    start: Option<isize>,
    /// An optional codon end coordinate in the gene to apply the rule to (1-based inclusive).
    /// If not given, the end of the gene is inferred
    end: Option<isize>,
    /// The drug(s) the rule applies to. Multiple drugs can be listed using a semi-colon delimiter
    #[serde(deserialize_with = "str_to_set")]
    drugs: BTreeSet<String>,
}

/// Allow serde to deserialize the drugs section of the expert rules
fn str_to_set<'de, D>(deserializer: D) -> Result<BTreeSet<String>, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    Ok(s.split(';').map(|item| item.to_owned()).collect())
}

#[cfg(test)]
mod tests {
    use crate::expert::{ExpertError, Rule, VariantType};
    use std::collections::BTreeSet;
    use std::str::FromStr;

    #[test]
    fn variant_type_display() {
        assert_eq!(VariantType::Frameshift.to_string(), "frameshift");
        assert_eq!(VariantType::Missense.to_string(), "missense");
        assert_eq!(VariantType::Nonsense.to_string(), "nonsense");
    }

    #[test]
    fn variant_type_from_str() {
        assert_eq!(VariantType::from_str("missense"), Ok(VariantType::Missense));
        assert_eq!(VariantType::from_str("nonsense"), Ok(VariantType::Nonsense));
        assert_eq!(
            VariantType::from_str("frameshift"),
            Ok(VariantType::Frameshift)
        );
        assert_eq!(
            VariantType::from_str("foo"),
            Err(ExpertError::UnknownVariantType("foo".to_string()))
        );
    }

    #[test]
    fn deserialise_correct_rule() {
        const CSV: &[u8] = b"frameshift,gene,1,10,drug1";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Rule = reader.deserialize().next().unwrap().unwrap();

        let expected = Rule {
            variant_type: VariantType::Frameshift,
            gene: "gene".to_string(),
            start: Some(1),
            end: Some(10),
            drugs: BTreeSet::from(["drug1".to_string()]),
        };

        assert_eq!(record, expected)
    }

    #[test]
    fn deserialise_correct_rule_multiple_drugs_sorted() {
        const CSV: &[u8] = b"frameshift,gene,1,10,drug1;drug2";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Rule = reader.deserialize().next().unwrap().unwrap();

        let expected = Rule {
            variant_type: VariantType::Frameshift,
            gene: "gene".to_string(),
            start: Some(1),
            end: Some(10),
            drugs: BTreeSet::from(["drug1".to_string(), "drug2".to_string()]),
        };

        assert_eq!(record, expected)
    }

    #[test]
    fn deserialise_correct_rule_multiple_drugs_unsorted() {
        const CSV: &[u8] = b"frameshift,gene,1,10,drug2;drug1";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Rule = reader.deserialize().next().unwrap().unwrap();

        let expected = Rule {
            variant_type: VariantType::Frameshift,
            gene: "gene".to_string(),
            start: Some(1),
            end: Some(10),
            drugs: BTreeSet::from(["drug1".to_string(), "drug2".to_string()]),
        };

        assert_eq!(record, expected)
    }

    #[test]
    fn deserialise_correct_rule_no_start_or_end() {
        const CSV: &[u8] = b"missense,gene,,,drug2";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Rule = reader.deserialize().next().unwrap().unwrap();

        let expected = Rule {
            variant_type: VariantType::Missense,
            gene: "gene".to_string(),
            start: None,
            end: None,
            drugs: BTreeSet::from(["drug2".to_string()]),
        };

        assert_eq!(record, expected)
    }

    #[test]
    fn deserialise_correct_rule_no_start() {
        const CSV: &[u8] = b"missense,gene,,10,drug2";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Rule = reader.deserialize().next().unwrap().unwrap();

        let expected = Rule {
            variant_type: VariantType::Missense,
            gene: "gene".to_string(),
            start: None,
            end: Some(10),
            drugs: BTreeSet::from(["drug2".to_string()]),
        };

        assert_eq!(record, expected)
    }

    #[test]
    fn deserialise_correct_rule_no_end() {
        const CSV: &[u8] = b"missense,gene,4,,drug2";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Rule = reader.deserialize().next().unwrap().unwrap();

        let expected = Rule {
            variant_type: VariantType::Missense,
            gene: "gene".to_string(),
            start: Some(4),
            end: None,
            drugs: BTreeSet::from(["drug2".to_string()]),
        };

        assert_eq!(record, expected)
    }

    #[test]
    fn deserialise_rule_wrong_delimiter() {
        const CSV: &[u8] = b"missense\tgene\t4\t5\tdrug2";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Result<Rule, csv::Error> = reader.deserialize().next().unwrap();
        assert!(record.is_err());
    }

    #[test]
    fn deserialise_rule_unknown_variant_type() {
        const CSV: &[u8] = b"foo,gene,4,5,drug2";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Result<Rule, csv::Error> = reader.deserialize().next().unwrap();
        assert!(record.is_err());
    }

    #[test]
    fn deserialise_rule_char_passed_for_start() {
        const CSV: &[u8] = b"missense,gene,s,5,drug2";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Result<Rule, csv::Error> = reader.deserialize().next().unwrap();
        assert!(record.is_err());
    }

    #[test]
    fn deserialise_rule_with_missing_gene() {
        const CSV: &[u8] = b"missense,4,5,drug2";

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(CSV);

        let record: Result<Rule, csv::Error> = reader.deserialize().next().unwrap();
        assert!(record.is_err());
    }
}
