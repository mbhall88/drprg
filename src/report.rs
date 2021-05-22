use crate::panel::{Residue, Variant};
use crate::predict::Prediction;
use serde::{Deserialize, Serialize};

#[derive(Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct Susceptibility {
    predict: Prediction,
    drug: String,
    evidence: Vec<Evidence>,
}
#[derive(Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct Evidence {
    variant: Variant,
    gene: String,
    residue: Residue,
    vcfid: String,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    fn remove_whitespace(s: &str) -> String {
        s.split_whitespace().collect()
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
            drug: "Ofloxacin".to_string(),
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
          "drug": "Ofloxacin",
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
