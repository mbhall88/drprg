use crate::report::Evidence;
use log::warn;
use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
use serde_derive::Deserialize;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fmt;
use std::hash::Hash;
use std::path::Path;
use std::str::FromStr;
use thiserror::Error;

#[derive(Error, Debug, PartialEq, Eq)]
pub enum ExpertError {
    /// Unrecognised variant type
    #[error("{0} is not a recognised variant type")]
    UnknownVariantType(String),
}

#[derive(Eq, PartialEq, Debug, Hash, Copy, Clone)]
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
#[derive(Debug, Deserialize, Eq, PartialEq, Hash, Clone)]
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
    pub(crate) drugs: BTreeSet<String>,
}

impl Rule {
    fn contains(&self, mutation: &Evidence) -> bool {
        if self.gene != mutation.gene {
            return false;
        }
        let start = self.start.unwrap_or(1);
        let end = self.end.unwrap_or(isize::MAX);
        let rule_range = start..=end;
        let mutation_pos = mutation.variant.pos as isize;
        if !rule_range.contains(&mutation_pos) {
            return false;
        }
        match self.variant_type {
            VariantType::Frameshift if mutation.is_frameshift() => true,
            VariantType::Missense if mutation.is_missense() => true,
            VariantType::Nonsense if mutation.is_nonsense() => true,
            _ => false,
        }
    }
}

/// Allow serde to deserialize the drugs section of the expert rules
fn str_to_set<'de, D>(deserializer: D) -> Result<BTreeSet<String>, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    Ok(s.split(';').map(|item| item.to_owned()).collect())
}

pub(crate) type ExpertRules = HashMap<String, HashSet<Rule>>;

pub trait RuleExt {
    fn from_csv(path: &Path) -> Result<ExpertRules, anyhow::Error>;
    fn matches(&self, mutation: &Evidence) -> Vec<Rule>;
}

impl RuleExt for ExpertRules {
    fn from_csv(path: &Path) -> Result<Self, anyhow::Error> {
        let mut rules: ExpertRules = HashMap::new();
        let mut n_records = 0;

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_path(path)?;

        for result in reader.deserialize() {
            n_records += 1;
            let record: Rule = result?;
            let set_of_records = rules
                .entry(record.gene.to_owned())
                .or_insert_with(HashSet::new);
            let seen_before = !set_of_records.insert(record);

            if seen_before {
                warn!(
                    "Duplicate expert rule detected at line number {}",
                    n_records
                )
            }
        }
        Ok(rules)
    }

    fn matches(&self, mutation: &Evidence) -> Vec<Rule> {
        let mut matches = Vec::new();
        if let Some(rules) = self.get(&mutation.gene) {
            for rule in rules {
                if rule.contains(mutation) {
                    matches.push(rule.to_owned());
                }
            }
        }
        matches
    }
}

#[cfg(test)]
mod tests {
    use crate::expert::{ExpertError, ExpertRules, Rule, RuleExt, VariantType};
    use crate::panel::{Residue, Variant};
    use crate::report::Evidence;
    use std::collections::{BTreeSet, HashMap, HashSet};
    use std::fs::File;
    use std::io::Write;
    use std::str::FromStr;
    use tempfile::NamedTempFile;

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

    #[test]
    fn expert_rules_from_csv() {
        const CSV: &[u8] = b"missense,geneA,1,2,drug\nnonsense,geneB,1,,drug;foo";
        let tmp = NamedTempFile::new().unwrap();
        {
            let mut file = File::create(tmp.path()).unwrap();
            file.write_all(CSV).unwrap();
        }
        let rules = ExpertRules::from_csv(tmp.path()).unwrap();
        let mut expected: ExpertRules = HashMap::new();
        expected.insert(
            "geneA".to_string(),
            HashSet::from([Rule {
                variant_type: VariantType::Missense,
                gene: "geneA".to_string(),
                start: Some(1),
                end: Some(2),
                drugs: BTreeSet::from(["drug".to_string()]),
            }]),
        );
        expected.insert(
            "geneB".to_string(),
            HashSet::from([Rule {
                variant_type: VariantType::Nonsense,
                gene: "geneB".to_string(),
                start: Some(1),
                end: None,
                drugs: BTreeSet::from(["foo".to_string(), "drug".to_string()]),
            }]),
        );

        assert_eq!(rules, expected)
    }

    #[test]
    fn expert_rules_from_csv_duplicate_rule_ignored() {
        const CSV: &[u8] = b"missense,geneA,1,2,drug\nmissense,geneA,1,2,drug";
        let tmp = NamedTempFile::new().unwrap();
        {
            let mut file = File::create(tmp.path()).unwrap();
            file.write_all(CSV).unwrap();
        }
        let rules = ExpertRules::from_csv(tmp.path()).unwrap();
        let mut expected: ExpertRules = HashMap::new();
        expected.insert(
            "geneA".to_string(),
            HashSet::from([Rule {
                variant_type: VariantType::Missense,
                gene: "geneA".to_string(),
                start: Some(1),
                end: Some(2),
                drugs: BTreeSet::from(["drug".to_string()]),
            }]),
        );

        assert_eq!(rules, expected)
    }

    #[test]
    fn expert_rules_from_csv_malformatted() {
        const CSV: &[u8] = b"missense,geneA,1,2,drug\nnonsense,geneB,1,drug;foo";
        let tmp = NamedTempFile::new().unwrap();
        {
            let mut file = File::create(tmp.path()).unwrap();
            file.write_all(CSV).unwrap();
        }
        let rules = ExpertRules::from_csv(tmp.path());
        assert!(rules.is_err());
    }

    #[test]
    fn rule_contains_gene_doesnt_match() {
        let mutation = Evidence {
            variant: Variant::from_str("A1T").unwrap(),
            gene: "foo".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Frameshift,
            gene: "bar".to_string(),
            start: None,
            end: None,
            drugs: Default::default(),
        };

        assert!(!rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_variant_types_differ() {
        let mutation = Evidence {
            variant: Variant::from_str("A1T").unwrap(),
            gene: "foo".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Frameshift,
            gene: "foo".to_string(),
            start: None,
            end: None,
            drugs: Default::default(),
        };

        assert!(!rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_upstream_of_start() {
        let mutation = Evidence {
            variant: Variant::from_str("A-1T").unwrap(),
            gene: "foo".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Missense,
            gene: "foo".to_string(),
            start: None,
            end: None,
            drugs: Default::default(),
        };

        assert!(!rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_at_start() {
        let mutation = Evidence {
            variant: Variant::from_str("A1T").unwrap(),
            gene: "foo".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Missense,
            gene: "foo".to_string(),
            start: None,
            end: None,
            drugs: Default::default(),
        };

        assert!(rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_past_end() {
        let mutation = Evidence {
            variant: Variant::from_str("A5T").unwrap(),
            gene: "foo".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Missense,
            gene: "foo".to_string(),
            start: None,
            end: Some(4),
            drugs: Default::default(),
        };

        assert!(!rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_at_end() {
        let mutation = Evidence {
            variant: Variant::from_str("A5T").unwrap(),
            gene: "foo".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Missense,
            gene: "foo".to_string(),
            start: None,
            end: Some(5),
            drugs: Default::default(),
        };

        assert!(rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_single_position() {
        let mutation = Evidence {
            variant: Variant::from_str("A5T").unwrap(),
            gene: "foo".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Missense,
            gene: "foo".to_string(),
            start: Some(5),
            end: Some(5),
            drugs: Default::default(),
        };

        assert!(rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_non_frameshift_indel() {
        let mutation = Evidence {
            variant: Variant::from_str("ACGT5A").unwrap(),
            gene: "foo".to_string(),
            residue: Default::default(),
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Frameshift,
            gene: "foo".to_string(),
            start: None,
            end: None,
            drugs: Default::default(),
        };

        assert!(!rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_frameshift_indel() {
        let mutation = Evidence {
            variant: Variant::from_str("ACG5A").unwrap(),
            gene: "foo".to_string(),
            residue: Default::default(),
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Frameshift,
            gene: "foo".to_string(),
            start: None,
            end: None,
            drugs: Default::default(),
        };

        assert!(rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_frameshift_in_promoter() {
        let mutation = Evidence {
            variant: Variant::from_str("ACG-5A").unwrap(),
            gene: "foo".to_string(),
            residue: Default::default(),
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Frameshift,
            gene: "foo".to_string(),
            start: Some(-7),
            end: None,
            drugs: Default::default(),
        };

        assert!(rule.contains(&mutation))
    }

    #[test]
    fn rule_contains_nonsense() {
        let mutation = Evidence {
            variant: Variant::from_str("A5*").unwrap(),
            gene: "foo".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };
        let rule = Rule {
            variant_type: VariantType::Nonsense,
            gene: "foo".to_string(),
            start: Some(5),
            end: None,
            drugs: Default::default(),
        };

        assert!(rule.contains(&mutation))
    }

    #[test]
    fn expert_rules_matches_variant_type_matches_but_wrong_gene() {
        const CSV: &[u8] = b"missense,geneA,1,2,drug\nnonsense,geneB,,,drug;foo";
        let tmp = NamedTempFile::new().unwrap();
        {
            let mut file = File::create(tmp.path()).unwrap();
            file.write_all(CSV).unwrap();
        }
        let rules = ExpertRules::from_csv(tmp.path()).unwrap();
        let mutation = Evidence {
            variant: Variant::from_str("A1T").unwrap(),
            gene: "geneB".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };

        assert!(rules.matches(&mutation).is_empty())
    }

    #[test]
    fn expert_rules_matches_vartype_and_gene_matches_but_not_range() {
        const CSV: &[u8] = b"missense,geneA,1,2,drug\nnonsense,geneB,,,drug;foo";
        let tmp = NamedTempFile::new().unwrap();
        {
            let mut file = File::create(tmp.path()).unwrap();
            file.write_all(CSV).unwrap();
        }
        let rules = ExpertRules::from_csv(tmp.path()).unwrap();
        let mutation = Evidence {
            variant: Variant::from_str("A3T").unwrap(),
            gene: "geneA".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };

        assert!(rules.matches(&mutation).is_empty())
    }

    #[test]
    fn expert_rules_matches() {
        const CSV: &[u8] = b"missense,geneA,1,2,drug\nnonsense,geneB,,,drug;foo";
        let tmp = NamedTempFile::new().unwrap();
        {
            let mut file = File::create(tmp.path()).unwrap();
            file.write_all(CSV).unwrap();
        }
        let rules = ExpertRules::from_csv(tmp.path()).unwrap();
        let mutation = Evidence {
            variant: Variant::from_str("A3*").unwrap(),
            gene: "geneB".to_string(),
            residue: Residue::Amino,
            vcfid: "".to_string(),
        };

        let actual = rules.matches(&mutation);
        let expected = vec![Rule {
            variant_type: VariantType::Nonsense,
            gene: "geneB".to_string(),
            start: None,
            end: None,
            drugs: BTreeSet::from(["drug".to_string(), "foo".to_string()]),
        }];

        assert_eq!(actual, expected)
    }
}
