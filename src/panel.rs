use std::collections::{HashMap, HashSet};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::str::FromStr;

use regex::Regex;
use rust_htslib::bcf;
use serde::{de, Deserialize, Deserializer};
use thiserror::Error;

pub(crate) type Panel = HashMap<String, HashSet<PanelRecord>>;
lazy_static! {
    static ref VARIANT_REGEX: Regex =
        Regex::new(r"^([a-zA-Z]+)(-?\d+)([a-zA-Z]+)$").unwrap();
}

/// A collection of custom errors relating to the command line interface for this package.
#[derive(Error, Debug, PartialEq)]
pub enum PanelError {
    /// Unrecognised residue type
    #[error("{0} is not a known residue type")]
    UnknownResidueType(String),
    /// The variant string is not in the correct format
    #[error("The variant is not in the correct format [<STR><INT><STR>]: {0}")]
    IncorrectVariantFormat(String),
    /// The record's gene is not in the VCF header
    #[error("Gene {0} is not in the VCF header")]
    GeneNotInHeader(String),
    /// Variant position is out of range for the gene and padding used
    #[error("The variant position {0} is out of range based on the padding and gene start for {1}")]
    PosOutOfRange(i64, String),
    /// Failed to set a VCF field
    #[error("Failed to set the VCF field {0} for {1}")]
    SetVcfFieldFailed(String, String),
    /// Reference allele of the panel record does not match the reference sequence
    #[error("Reference allele for {0} does not match the reference sequence")]
    RefDoesNotMatch(String),
    /// Got more than one amino acid in variant
    #[error("No support for multiple amino acid allele variants [{0}]")]
    MultiAminoAlleleNotSupported(String),
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

impl PanelRecord {
    /// The unique name for this record, which is the joining of the gene and variant
    pub fn name(&self) -> String {
        format!("{}_{}", self.gene, self.variant)
    }
    /// The reference allele for the variant this record describes
    fn ref_allele(&self) -> &[u8] {
        self.variant.reference.as_bytes()
    }
    /// Returns all possible reference alleles; converting amino acids to all codons if necessary
    fn all_ref_alleles(&self) -> Result<Vec<&[u8]>, PanelError> {
        match self.residue {
            Residue::Nucleic => Ok(vec![self.ref_allele()]),
            Residue::Amino => {
                if self.ref_allele().len() > 1 {
                    Err(PanelError::MultiAminoAlleleNotSupported(self.name()))
                } else {
                    Ok(amino_to_codons(self.ref_allele()[0]))
                }
            }
        }
    }
    /// The position of the variant within the gene/protein
    fn pos(&self) -> i64 {
        self.variant.pos
    }
    /// The alternate allele the variant describes
    fn alt_allele(&self) -> &[u8] {
        self.variant.new.as_bytes()
    }
    /// Generate the header entries for the panel record INFO fields
    pub fn vcf_header_entries() -> Vec<&'static [u8]> {
        vec![
            b"##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene the variant occurs in\">",
            b"##INFO=<ID=VAR,Number=1,Type=String,Description=\"The variant describing reference, position, alternate on the gene\">",
            b"##INFO=<ID=RES,Number=1,Type=String,Description=\"Residue the variant describes (i.e. Nucleic/Amino)\">",
            b"##INFO=<ID=DRUGS,Number=.,Type=String,Description=\"Drugs this variant causes resistance to\">",
            br#"##INFO=<ID=PAD,Number=1,Type=Integer,Description="Number of bases added to start and end of gene">"#,
            br#"##INFO=<ID=ST,Number=1,Type=String,Description="Strand the gene is on">"#,
        ]
    }
    /// Check the panel reference allele matches the reference sequence, If it does, output the DNA
    /// version of the reference.
    ///
    /// # Errors
    /// Returns a [`PanelError::RefDoesNotMatch`] if no reference allele matches the reference sequence.
    fn check_ref(&self, refseq: &Vec<u8>, padding: u32) -> Result<&[u8], PanelError> {
        let pos = self.pos() + padding as i64;
        let ref_alleles = self.all_ref_alleles();
        Ok(b"dummy")
    }

    /// Fill a VCF `record` with information about this panel record.
    pub fn to_vcf(
        &self,
        record: &mut bcf::Record,
        padding: u32,
    ) -> Result<(), PanelError> {
        let pos = self.pos() + (padding as i64) - 1; // rust htslib works with 0-based pos
        if pos < 0 {
            return Err(PanelError::PosOutOfRange(self.pos(), self.gene.to_owned()));
        }
        let rid = match record.header().name2rid(self.gene.as_bytes()) {
            Ok(i) => i,
            Err(_) => return Err(PanelError::GeneNotInHeader(self.gene.to_owned())),
        };
        record.set_rid(Some(rid));
        record.set_pos(pos);
        record
            .set_id(self.name().as_bytes())
            .map_err(|_| PanelError::SetVcfFieldFailed("ID".to_owned(), self.name()))?;
        record
            .push_info_integer(b"PAD", [padding as i32].as_ref())
            .map_err(|_| {
                PanelError::SetVcfFieldFailed("PAD".to_owned(), self.name())
            })?;
        record
            .push_info_string(b"GENE", &[self.gene.as_bytes()])
            .map_err(|_| {
                PanelError::SetVcfFieldFailed("GENE".to_owned(), self.name())
            })?;
        record
            .push_info_string(b"VAR", &[self.variant.to_string().as_bytes()])
            .map_err(|_| {
                PanelError::SetVcfFieldFailed("VAR".to_owned(), self.name())
            })?;
        record
            .push_info_string(b"RES", &[self.residue.to_string().as_bytes()])
            .map_err(|_| {
                PanelError::SetVcfFieldFailed("RES".to_owned(), self.name())
            })?;
        let mut drugs: Vec<&[u8]> = self.drugs.iter().map(|d| d.as_bytes()).collect();
        drugs.sort();
        record
            .push_info_string(b"DRUGS", drugs.as_slice())
            .map_err(|_| {
                PanelError::SetVcfFieldFailed("DRUGS".to_owned(), self.name())
            })?;
        // todo: check and set ref
        // todo: set alts
        Ok(())
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

fn amino_to_codons(amino_acid: u8) -> Vec<&'static [u8]> {
    match amino_acid {
        b'F' => vec![b"TTT", b"TTC"], // (Phe/F) Phenylalanine
        b'L' => vec![b"TTA", b"TTG", b"CTT", b"CTC", b"CTA", b"CTG"], // (Leu/L) Leucine
        b'I' => vec![b"ATT", b"ATC", b"ATA"], // (Ile/I) Isoleucine
        b'M' => vec![b"ATG"],         // (Met/M) Methionine
        b'V' => vec![b"GTT", b"GTC", b"GTA", b"GTG"], // (Val/V) Valine
        b'S' => vec![b"TCT", b"TCC", b"TCA", b"TCG", b"AGT", b"AGC"], // (Ser/S) Serine
        b'P' => vec![b"CCT", b"CCC", b"CCA", b"CCG"], // (Pro/P) Proline
        b'T' => vec![b"ACT", b"ACC", b"ACA", b"ACG"], // (Thr/T) Threonine
        b'A' => vec![b"GCT", b"GCC", b"GCA", b"GCG"], // (Ala/A) Alanine
        b'Y' => vec![b"TAT", b"TAC"], // (Tyr/Y) Tyrosine
        b'H' => vec![b"CAT", b"CAC"], // (His/H) Histidine
        b'Q' => vec![b"CAA", b"CAG"], // (Gln/Q) Glutamine
        b'N' => vec![b"AAT", b"AAC"], // (Asn/N) Asparagine
        b'K' => vec![b"AAA", b"AAG"], // (Lys/K) Lysine
        b'D' => vec![b"GAT", b"GAC"], // (Asp/D) Aspartic acid
        b'E' => vec![b"GAA", b"GAG"], // (Glu/E) Glutamic acid
        b'C' => vec![b"TGT", b"TGC"], // (Cys/C) Cysteine
        b'W' => vec![b"TGG"],         // (Trp/W) Tryptophan
        b'R' => vec![b"CGT", b"CGC", b"CGA", b"CGG", b"AGA", b"AGG"], // (Arg/R) Arginine
        b'G' => vec![b"GGT", b"GGC", b"GGA", b"GGG"],                 // (Gly/G) Glycine
        b'*' => vec![b"TGA", b"TAA", b"TAG"], // Stop/Termination
        _ => vec![],
    }
}

#[cfg(test)]
mod tests {
    use std::iter::FromIterator;

    use super::*;
    use rust_htslib::bcf::{Format, Header, Writer};
    use std::ops::Deref;
    use tempfile::NamedTempFile;

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

    #[test]
    fn panel_record_vcf_header_entries() {
        let actual = PanelRecord::vcf_header_entries();
        let expected = vec![
            "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene the variant occurs in\">".as_bytes(),
            "##INFO=<ID=VAR,Number=1,Type=String,Description=\"The variant describing reference, position, alternate on the gene\">".as_bytes(),
            "##INFO=<ID=RES,Number=1,Type=String,Description=\"Residue the variant describes (i.e. Nucleic/Amino)\">".as_bytes(),
            "##INFO=<ID=DRUGS,Number=.,Type=String,Description=\"Drugs this variant causes resistance to\">".as_bytes(),
            br#"##INFO=<ID=PAD,Number=1,Type=Integer,Description="Number of bases added to start and end of gene">"#,
            br#"##INFO=<ID=ST,Number=1,Type=String,Description="Strand the gene is on">"#,
        ];
        assert_eq!(actual, expected)
    }

    #[test]
    fn panel_record_name_getter() {
        let record = PanelRecord {
            gene: "gene".to_string(),
            variant: Variant::from_str("K1S").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        assert_eq!(record.name(), "gene_K1S")
    }

    #[test]
    fn panel_record_ref_allele_getter() {
        let record = PanelRecord {
            gene: "gene".to_string(),
            variant: Variant::from_str("K1S").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        assert_eq!(record.ref_allele(), b"K")
    }

    #[test]
    fn panel_record_alt_allele_getter() {
        let record = PanelRecord {
            gene: "gene".to_string(),
            variant: Variant::from_str("K1S").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        assert_eq!(record.alt_allele(), b"S")
    }

    #[test]
    fn panel_record_pos_getter() {
        let record = PanelRecord {
            gene: "gene".to_string(),
            variant: Variant::from_str("K1S").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        assert_eq!(record.pos(), 1)
    }

    #[test]
    fn panel_record_to_vcf_pos_out_of_range() {
        let record = PanelRecord {
            gene: "gene".to_string(),
            variant: Variant::from_str("A-1T").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };
        let padding = 0;
        let mut header = Header::new();
        header.push_record(br#"##contig=<ID=gene,length=10>"#);
        let tmpfile = NamedTempFile::new().unwrap();
        let tmppath = tmpfile.path();
        let vcf = Writer::from_path(tmppath, &header, true, Format::VCF).unwrap();
        let mut vcf_record = vcf.empty_record();

        let actual = record.to_vcf(&mut vcf_record, padding).unwrap_err();
        let expected = PanelError::PosOutOfRange(-1, "gene".to_string());
        assert_eq!(actual, expected)
    }

    #[test]
    fn panel_record_to_vcf_gene_not_in_header() {
        let record = PanelRecord {
            gene: "gene".to_string(),
            variant: Variant::from_str("A1T").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };
        let padding = 0;
        let mut header = Header::new();
        header.push_record(r#"##contig=<ID=foo,length=10>"#.as_bytes());
        let tmpfile = NamedTempFile::new().unwrap();
        let tmppath = tmpfile.path();
        let vcf = Writer::from_path(tmppath, &header, true, Format::VCF).unwrap();
        let mut vcf_record = vcf.empty_record();

        let actual = record.to_vcf(&mut vcf_record, padding).unwrap_err();
        let expected = PanelError::GeneNotInHeader("gene".to_string());
        assert_eq!(actual, expected)
    }

    #[test]
    fn panel_record_to_vcf_with_two_drugs() {
        let record = PanelRecord {
            gene: "gene".to_string(),
            variant: Variant::from_str("A1T").unwrap(),
            residue: Residue::Nucleic,
            drugs: HashSet::from_iter(vec!["d1".to_string(), "d2".to_string()]),
        };
        let padding = 0;
        let mut header = Header::new();
        header.push_record(r#"##contig=<ID=gene,length=10>"#.as_bytes());
        for entry in PanelRecord::vcf_header_entries() {
            header.push_record(entry);
        }
        let tmpfile = NamedTempFile::new().unwrap();
        let tmppath = tmpfile.path();
        let vcf = Writer::from_path(tmppath, &header, true, Format::VCF).unwrap();
        let mut vcf_record = vcf.empty_record();

        record.to_vcf(&mut vcf_record, padding).unwrap();
        let actual = vcf_record.info(b"DRUGS").string().unwrap().unwrap();
        let expected = &[b"d1", b"d2"];
        assert_eq!(actual.deref(), expected)
    }

    #[test]
    fn amino_to_codons_serine() {
        let actual = amino_to_codons(b'S');
        let expected = vec![b"TCT", b"TCC", b"TCA", b"TCG", b"AGT", b"AGC"];

        assert_eq!(actual, expected)
    }

    #[test]
    fn amino_to_codons_stop() {
        let actual = amino_to_codons(b'*');
        let expected = vec![b"TGA", b"TAA", b"TAG"];

        assert_eq!(actual, expected)
    }

    #[test]
    fn amino_to_codons_unknown() {
        let actual = amino_to_codons(b'Z');

        assert!(actual.is_empty())
    }

    #[test]
    fn all_ref_alleles_nucleic_returns_self() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("CC3A").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };

        let actual = record.all_ref_alleles().unwrap();
        let expected = vec![b"CC"];
        assert_eq!(actual, expected)
    }

    #[test]
    fn all_ref_alleles_amino_returns_codons() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("C3A").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        let actual = record.all_ref_alleles().unwrap();
        let expected = vec![b"TGT", b"TGC"];
        assert_eq!(actual, expected)
    }

    #[test]
    fn all_ref_alleles_multi_amino_returns_err() {
        let record = PanelRecord {
            gene: "G".to_string(),
            variant: Variant::from_str("CW3A").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        let actual = record.all_ref_alleles().unwrap_err();
        let expected = PanelError::MultiAminoAlleleNotSupported("G_CW3A".to_string());
        assert_eq!(actual, expected)
    }
}
