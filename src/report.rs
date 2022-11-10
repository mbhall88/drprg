use crate::panel::{Residue, Variant};
use crate::predict::Prediction;
use serde_derive::{Deserialize, Serialize};
use std::str::FromStr;

const STOP: &str = "*";

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
    pub fn to_variant_string(&self) -> String {
        format!("{}_{}", self.gene, self.variant)
    }
    pub fn is_synonymous(&self) -> bool {
        self.residue == Residue::Amino && self.variant.reference == self.variant.new
    }
    pub fn is_missense(&self) -> bool {
        self.residue == Residue::Amino && !self.is_nonsense() && !self.is_synonymous()
    }
    pub fn is_nonsense(&self) -> bool {
        self.variant.new == STOP && self.residue == Residue::Amino
    }
    pub fn is_frameshift(&self) -> bool {
        let len_diff = self
            .variant
            .reference
            .len()
            .abs_diff(self.variant.new.len());
        self.residue == Residue::Nucleic && len_diff % 3 != 0
    }
    pub fn atomise(&self) -> Vec<Self> {
        if self.variant.is_snp() || self.variant.is_indel() {
            return vec![self.to_owned()];
        }

        let mut v = vec![];

        for (i, (r, a)) in self
            .variant
            .reference
            .chars()
            .zip(self.variant.new.chars())
            .enumerate()
        {
            // we unwrap here as we know we have valid data
            let var = Variant::from_str(
                format!("{}{}{}", r, self.variant.pos + i as i64, a).as_str(),
            )
            .unwrap();
            let ev = Evidence {
                variant: var,
                gene: self.gene.to_owned(),
                residue: self.residue.to_owned(),
                vcfid: self.vcfid.to_owned(),
            };
            v.push(ev);
        }

        v
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
    fn evidence_atomise_snp_returns_vec_of_same() {
        let ev = Evidence {
            variant: Variant::from_str("A4A").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };

        let actual = ev.atomise();
        let expected = vec![ev];

        assert_eq!(actual, expected)
    }
    #[test]
    fn evidence_atomise_mnp_returns_vec_of_snps() {
        let ev = Evidence {
            variant: Variant::from_str("AG4AT").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };

        let actual = ev.atomise();
        let expected = vec![
            Evidence {
                variant: Variant::from_str("A4A").unwrap(),
                gene: "inhA".to_string(),
                residue: Residue::Nucleic,
                vcfid: "abcd1234".to_string(),
            },
            Evidence {
                variant: Variant::from_str("G5T").unwrap(),
                gene: "inhA".to_string(),
                residue: Residue::Nucleic,
                vcfid: "abcd1234".to_string(),
            },
        ];

        assert_eq!(actual, expected)
    }
    #[test]
    fn evidence_atomise_indel_returns_vec_of_same() {
        let ev = Evidence {
            variant: Variant::from_str("A4CA").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };

        let actual = ev.atomise();
        let expected = vec![ev];

        assert_eq!(actual, expected)
    }
    #[test]
    fn evidence_atomise_single_amino_change_returns_vec_of_same() {
        let ev = Evidence {
            variant: Variant::from_str("D94G").unwrap(),
            gene: "gyrA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };

        let actual = ev.atomise();
        let expected = vec![ev];

        assert_eq!(actual, expected)
    }
    #[test]
    fn evidence_atomise_multi_amino_change_returns_vec_of_single_amino_changes() {
        let ev = Evidence {
            variant: Variant::from_str("DS94GT").unwrap(),
            gene: "gyrA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };

        let actual = ev.atomise();
        let expected = vec![
            Evidence {
                variant: Variant::from_str("D94G").unwrap(),
                gene: "gyrA".to_string(),
                residue: Residue::Amino,
                vcfid: "abcd1234".to_string(),
            },
            Evidence {
                variant: Variant::from_str("S95T").unwrap(),
                gene: "gyrA".to_string(),
                residue: Residue::Amino,
                vcfid: "abcd1234".to_string(),
            },
        ];

        assert_eq!(actual, expected)
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
    fn evidence_is_missense_nucleic_acid() {
        let ev = Evidence {
            variant: Variant::from_str("A4T").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_missense())
    }

    #[test]
    fn evidence_is_missense() {
        let ev = Evidence {
            variant: Variant::from_str("A4T").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        assert!(ev.is_missense())
    }

    #[test]
    fn evidence_is_missense_is_nonsense() {
        let ev = Evidence {
            variant: Variant::from_str("A4*").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_missense())
    }

    #[test]
    fn evidence_is_nonsense() {
        let ev = Evidence {
            variant: Variant::from_str("A4*").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        assert!(ev.is_nonsense())
    }

    #[test]
    fn evidence_is_nonsense_is_nonsense() {
        let ev = Evidence {
            variant: Variant::from_str("A4T").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_nonsense())
    }

    #[test]
    fn evidence_is_nonsense_is_synonymous() {
        let ev = Evidence {
            variant: Variant::from_str("A4A").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_nonsense())
    }

    #[test]
    fn evidence_is_nonsense_is_nucleic() {
        let ev = Evidence {
            variant: Variant::from_str("A4*").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_nonsense())
    }

    #[test]
    fn evidence_is_frameshift_is_snp() {
        let ev = Evidence {
            variant: Variant::from_str("A4T").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_frameshift())
    }

    #[test]
    fn evidence_is_frameshift_is_1bp_frameshift() {
        let ev = Evidence {
            variant: Variant::from_str("AA4T").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };
        assert!(ev.is_frameshift())
    }

    #[test]
    fn evidence_is_frameshift_is_2bp_frameshift() {
        let ev = Evidence {
            variant: Variant::from_str("AAA4T").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };
        assert!(ev.is_frameshift())
    }

    #[test]
    fn evidence_is_frameshift_is_3bp_indel() {
        let ev = Evidence {
            variant: Variant::from_str("AAAA4T").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };
        assert!(!ev.is_frameshift())
    }

    #[test]
    fn evidence_is_frameshift_is_4bp_indel() {
        let ev = Evidence {
            variant: Variant::from_str("AAAA4TAAAAAAA").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Nucleic,
            vcfid: "abcd1234".to_string(),
        };
        assert!(ev.is_frameshift())
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

    #[test]
    fn evidence_to_variant_str() {
        let ev = Evidence {
            variant: Variant::from_str("K4S").unwrap(),
            gene: "inhA".to_string(),
            residue: Residue::Amino,
            vcfid: "abcd1234".to_string(),
        };

        let actual = ev.to_variant_string();
        let expected = String::from("inhA_K4S");

        assert_eq!(actual, expected)
    }
}
