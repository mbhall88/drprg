use crate::panel::{Residue, Variant};
use crate::predict::Prediction;
use serde_derive::{Deserialize, Serialize};

#[derive(Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct Susceptibility {
    pub(crate) predict: Prediction,
    pub(crate) evidence: Vec<Evidence>,
}
#[derive(Debug, Default, Serialize, Deserialize, PartialEq, Eq, Clone)]
pub struct Evidence {
    pub(crate) variant: Variant,
    pub(crate) gene: String,
    pub(crate) residue: Residue,
    pub(crate) vcfid: String,
}

impl Evidence {
    pub fn is_synonymous(&self) -> bool {
        self.residue == Residue::Amino && self.variant.reference == self.variant.new
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    fn remove_whitespace(s: &str) -> String {
        s.split_whitespace().collect()
    }

    #[test]
    fn evidence_is_synonymous_nucelic_is_not() {
        let ev = Evidence {
            variant: Variant::from_str("A4A").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_synonymous())
    }

    #[test]
    fn evidence_is_synonymous_amino_is_not() {
        let ev = Evidence {
            variant: Variant::from_str("A4T").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_synonymous())
    }

    #[test]
    fn evidence_is_synonymous_amino_is() {
        let ev = Evidence {
            variant: Variant::from_str("A4A").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        assert!(ev.is_synonymous())
    }

    #[test]
    fn evidence_is_synonymous_amino_is_multi_base() {
        let ev = Evidence {
            variant: Variant::from_str("AT4AT").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        assert!(ev.is_synonymous())
    }

    #[test]
    fn evidence_is_synonymous_amino_is_not_multi_base() {
        let ev = Evidence {
            variant: Variant::from_str("AT4AC").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_synonymous())
    }

    #[test]
    fn evidence_serde() {
        let expected_struct = Evidence {
            variant: Variant::from_str("K4S").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        let actual_data = serde_json::to_string(&expected_struct).unwrap();
        // Some JSON input data as a &str. Maybe this comes from the user.
        let expected_data = remove_whitespace(
            r#"
        {
            "variant": "K4S",
            "gene": "inhA",
            "residue": "PROT",
            "vcfid": "abcd1234"
        }"#,
        );
        let actual_struct: Evidence = serde_json::from_str(&expected_data).unwrap();

        assert_eq!(actual_data, expected_data);
        assert_eq!(actual_struct, expected_struct)
    }

    #[test]
    fn susceptibility_serde() {
        let expected_struct = Susceptibility {
            predict: Prediction::Resistant,
            evidence: vec![Evidence {
                variant: Variant::from_str("K4S").unwrap(),
                gene: "inhA".to_string(),
                residue: Residue::Amino,
                vcfid: "abcd1234".to_string(),
            }],
        };
        let actual_data = serde_json::to_string(&expected_struct).unwrap();
        // Some JSON input data as a &str. Maybe this comes from the user.
        let expected_data = remove_whitespace(
            r#"
        {
          "predict": "R",
          "evidence": [
            {
              "variant": "K4S",
              "gene": "inhA",
              "residue": "PROT",
              "vcfid": "abcd1234"
            }
          ]
        }"#,
        );
        let actual_struct: Susceptibility =
            serde_json::from_str(&expected_data).unwrap();

        assert_eq!("R", Prediction::Resistant.to_string());
        assert_eq!(actual_data, expected_data);
        assert_eq!(actual_struct, expected_struct)
    }
}
