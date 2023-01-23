use log::warn;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::ops::RangeInclusive;
use std::str::FromStr;

use crate::report::STOP;
use bstr::ByteSlice;
use regex::Regex;
use rust_htslib::bcf;
use serde::{self, de, Deserialize, Deserializer, Serialize, Serializer};
use serde_derive::Deserialize;
use std::path::Path;
use thiserror::Error;

pub(crate) type Panel = HashMap<String, HashSet<PanelRecord>>;

pub trait PanelExt {
    fn from_csv(path: &Path) -> Result<Panel, anyhow::Error>;
}

impl PanelExt for Panel {
    fn from_csv(path: &Path) -> Result<Self, anyhow::Error> {
        let mut panel: Panel = HashMap::new();
        let mut n_records = 0;

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(path)?;

        for result in reader.deserialize() {
            n_records += 1;
            let record: PanelRecord = result?;
            let set_of_records = panel
                .entry(record.gene.to_owned())
                .or_insert_with(HashSet::new);
            let seen_before = !set_of_records.insert(record);

            if seen_before {
                warn!(
                    "Duplicate panel record detected in record number {}",
                    n_records
                )
            }
        }
        Ok(panel)
    }
}

lazy_static! {
    static ref VARIANT_REGEX: Regex =
        Regex::new(r"^([a-zA-Z\*]+)(-?\d+)([a-zA-Z\*]+)$").unwrap();
    static ref NUCLEOTIDES: Vec<&'static [u8]> = vec![b"A", b"C", b"G", b"T"];
}
static AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

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
    /// Represents all other htslib errors
    #[error(transparent)]
    BcfError(#[from] rust_htslib::errors::Error),
    /// Got a negative position for an amino acid
    #[error("Negative positions are not allowed for protein residues [{0}]")]
    NegativeAminoPosition(String),
}

/// An enum representing the panel residue types we recognise
#[derive(Debug, PartialEq, Hash, Eq, Clone)]
pub enum Residue {
    Nucleic,
    Amino,
}

impl Default for Residue {
    fn default() -> Self {
        Self::Nucleic
    }
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

impl Serialize for Residue {
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
#[derive(Debug, Default, PartialEq, Hash, Eq, Clone)]
pub struct Variant {
    pub reference: String,
    pub pos: i64,
    pub new: String,
}

impl Variant {
    pub fn simplify(&self) -> Self {
        if self.reference != self.new {
            let mut reference = self.reference.to_owned();
            let mut new = self.new.to_owned();

            let mut pos = self.pos.to_owned();

            while reference.chars().next() == new.chars().next()
                && reference.len() != 1
                && new.len() != 1
            {
                reference.remove(0);
                new.remove(0);
                pos += 1;
            }
            while reference.chars().last() == new.chars().last()
                && reference.len() != 1
                && new.len() != 1
            {
                reference.pop();
                new.pop();
            }
            Variant {
                reference,
                pos,
                new,
            }
        } else {
            self.clone()
        }
    }
    pub fn is_indel(&self) -> bool {
        self.new.len() != self.reference.len()
    }
    pub fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.new.len() == 1
    }
    pub fn gene_deletion() -> Self {
        Self {
            reference: "".to_string(),
            pos: 0,
            new: "-".to_string(),
        }
    }
    pub fn start_lost() -> Self {
        Self {
            reference: "".to_string(),
            pos: 1,
            new: "-".to_string(),
        }
    }
    pub fn stop_lost(pos: i64) -> Self {
        Self {
            reference: STOP.to_string(),
            pos,
            new: "-".to_string(),
        }
    }
    pub fn is_gene_deletion(&self) -> bool {
        self.reference.is_empty() && self.pos == 0 && self.new == "-"
    }
    pub fn is_start_lost(&self) -> bool {
        self.reference.is_empty() && self.pos == 1 && self.new == "-"
    }
    pub fn is_stop_lost(&self) -> bool {
        self.reference == STOP && self.pos >= 1 && self.new == "-"
    }

    pub fn range(&self) -> RangeInclusive<i64> {
        let len = self.reference.len() as i64;
        let mut end = self.pos + (len - 1);
        if self.pos.is_negative() && end > -1 {
            end += 1;
        }
        self.pos..=end
    }
}

impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let to_write = if self.is_gene_deletion() {
            "gene_absent".to_string()
        } else if self.is_start_lost() {
            "start_lost".to_string()
        } else if self.is_stop_lost() {
            "stop_lost".to_string()
        } else {
            format!("{}{}{}", &self.reference, self.pos, &self.new)
        };
        write!(f, "{}", to_write)
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

impl Serialize for Variant {
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
    /// The position within the gene. That is, it converts protein position to DNA position.
    fn gene_pos(&self) -> Result<i64, PanelError> {
        let pos = match self.residue {
            Residue::Nucleic => {
                if self.pos() < 1 {
                    self.pos() + 1
                } else {
                    self.pos()
                }
            }
            Residue::Amino => {
                if self.pos() < 1 {
                    return Err(PanelError::NegativeAminoPosition(self.name()));
                } else {
                    3 * self.pos() - 2
                }
            }
        };
        Ok(pos)
    }
    /// The alternate allele the variant describes
    fn alt_allele(&self) -> &[u8] {
        self.variant.new.as_bytes()
    }
    /// Returns all possible alternate alleles; converting amino acids into codons if necessary.
    fn all_alt_alleles(&self) -> Result<Vec<&[u8]>, PanelError> {
        if !self.alt_allele().contains(&b'X') {
            match self.residue {
                Residue::Nucleic => Ok(vec![self.alt_allele()]),
                Residue::Amino => {
                    if self.alt_allele().len() > 1 {
                        Err(PanelError::MultiAminoAlleleNotSupported(self.name()))
                    } else {
                        Ok(amino_to_codons(self.alt_allele()[0]))
                    }
                }
            }
        } else if self.alt_allele().len() > 1 {
            Err(PanelError::MultiAminoAlleleNotSupported(self.name()))
        } else {
            match self.residue {
                Residue::Nucleic => Ok(NUCLEOTIDES
                    .iter()
                    .filter(|c| *c != &self.ref_allele())
                    .map(|c| c.to_owned())
                    .collect()),
                Residue::Amino => Ok(AMINO_ACIDS
                    .iter()
                    .filter(|&c| *c != self.ref_allele()[0])
                    .flat_map(|c| amino_to_codons(*c))
                    .collect()),
            }
        }
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
    fn check_ref<'a>(
        &self,
        refseq: &'a [u8],
        padding: u32,
    ) -> Result<&'a [u8], PanelError> {
        let ref_alleles = self.all_ref_alleles()?;
        if ref_alleles.is_empty() {
            return Err(PanelError::RefDoesNotMatch(self.name()));
        }
        let ref_len = ref_alleles[0].len();
        let start = (self.gene_pos()? - 1 + padding as i64) as usize;
        let end = start + ref_len;
        let expected_ref: &[u8] = &refseq[start..end];
        match ref_alleles
            .iter()
            .position(|&allele| allele == expected_ref)
        {
            Some(_) => Ok(expected_ref.as_bytes()),
            None => Err(PanelError::RefDoesNotMatch(self.name())),
        }
    }

    /// Fill a VCF `record` with information about this panel record.
    pub fn to_vcf(
        &self,
        record: &mut bcf::Record,
        refseq: &[u8],
        padding: u32,
    ) -> Result<(), PanelError> {
        let pos = self.gene_pos()? + (padding as i64) - 1; // rust htslib works with 0-based pos
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
        let ref_allele = self.check_ref(refseq, padding)?;
        let alt_alleles = self.all_alt_alleles()?;
        let mut alleles = vec![ref_allele];
        alleles.extend(alt_alleles);
        record.set_alleles(&alleles)?;
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
    use std::ops::Deref;

    use rust_htslib::bcf::{Format, Header, Writer};
    use tempfile::NamedTempFile;

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
    fn variant_from_str_with_stop_codon() {
        let s = "K2*";
        let actual = Variant::from_str(s).unwrap();
        let expected = Variant {
            reference: "K".to_string(),
            pos: 2,
            new: "*".to_string(),
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
    fn variant_range_snp() {
        let s = Variant::from_str("K2*").unwrap();
        let actual = s.range();
        let expected = 2..=2;

        assert_eq!(actual, expected)
    }

    #[test]
    fn variant_range_mnp() {
        let s = Variant::from_str("ATC2TTC").unwrap();
        let actual = s.range();
        let expected = 2..=4;

        assert_eq!(actual, expected)
    }

    #[test]
    fn variant_range_indel() {
        let s = Variant::from_str("ATC2TC").unwrap();
        let actual = s.range();
        let expected = 2..=4;

        assert_eq!(actual, expected)
    }

    #[test]
    fn variant_range_indel_in_promoter() {
        let s = Variant::from_str("ATC-6TC").unwrap();
        let actual = s.range();
        let expected = -6..=-4;

        assert_eq!(actual, expected)
    }

    #[test]
    fn variant_range_indel_in_promoter_crosses_start_pos() {
        let s = Variant::from_str("ATC-2TC").unwrap();
        let actual = s.range();
        let expected = -2..=1;

        assert_eq!(actual, expected)
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
    fn panel_record_gene_pos_for_nucleic_is_same() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("C6A").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };

        let actual = record.gene_pos().unwrap();
        let expected = 6;

        assert_eq!(actual, expected)
    }

    #[test]
    fn panel_record_gene_pos_for_amino_is_adjusted() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("C6A").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        let actual = record.gene_pos().unwrap();
        let expected = 16;

        assert_eq!(actual, expected)
    }

    #[test]
    fn panel_record_gene_pos_for_first_amino_is_one() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("C1A").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        let actual = record.gene_pos().unwrap();
        let expected = 1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn panel_record_gene_pos_for_negative_one_is_adjusted_to_zero() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("C-1A").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };

        let actual = record.gene_pos().unwrap();
        let expected = 0;

        assert_eq!(actual, expected)
    }

    #[test]
    fn panel_record_gene_pos_for_negative_is_adjusted() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("C-12A").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };

        let actual = record.gene_pos().unwrap();
        let expected = -11;

        assert_eq!(actual, expected)
    }

    #[test]
    fn panel_record_gene_pos_negative_amino_returns_err() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("C-12A").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        let actual = record.gene_pos().unwrap_err();
        let expected = PanelError::NegativeAminoPosition(record.name());

        assert_eq!(actual, expected)
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
        let vcf = Writer::from_path(tmppath, &header, true, Format::Vcf).unwrap();
        let mut vcf_record = vcf.empty_record();

        let actual = record.to_vcf(&mut vcf_record, &[], padding).unwrap_err();
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
        let vcf = Writer::from_path(tmppath, &header, true, Format::Vcf).unwrap();
        let mut vcf_record = vcf.empty_record();

        let actual = record.to_vcf(&mut vcf_record, &[], padding).unwrap_err();
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
        let vcf = Writer::from_path(tmppath, &header, true, Format::Vcf).unwrap();
        let mut vcf_record = vcf.empty_record();
        let refseq = b"A".to_vec();

        record.to_vcf(&mut vcf_record, &refseq, padding).unwrap();
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

    #[test]
    fn check_ref_nucleic_ref_matches() {
        let record = PanelRecord {
            gene: "G".to_string(),
            variant: Variant::from_str("CC1A").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };
        let refseq = b"AACCTTGG".to_vec();
        let padding = 2;

        let actual = record.check_ref(&refseq, padding).unwrap();
        let expected = b"CC";

        assert_eq!(actual, expected)
    }

    #[test]
    fn check_ref_nucleic_ref_does_not_match() {
        let record = PanelRecord {
            gene: "G".to_string(),
            variant: Variant::from_str("CC2A").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };
        let refseq = b"AACCTTGG".to_vec();
        let padding = 2;

        let actual = record.check_ref(&refseq, padding).unwrap_err();
        let expected = PanelError::RefDoesNotMatch(record.name());

        assert_eq!(actual, expected)
    }

    #[test]
    fn check_ref_nucleic_ref_negative_pos() {
        let record = PanelRecord {
            gene: "G".to_string(),
            variant: Variant::from_str("T-12C").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };
        let refseq = b"ACGTATGGTGGACGTATGCGGGCGTTGATC".to_vec();
        let padding = 15;

        let actual = record.check_ref(&refseq, padding).unwrap();
        let expected = b"T";

        assert_eq!(actual, expected)
    }

    #[test]
    fn check_ref_nucleic_ref_negative_pos_multiple_bases() {
        let record = PanelRecord {
            gene: "G".to_string(),
            variant: Variant::from_str("TTT-12C").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };
        let refseq = b"ACGTTTGGTGGACGTATGCGGGCGTTGATC".to_vec();
        let padding = 15;

        let actual = record.check_ref(&refseq, padding).unwrap();
        let expected = b"TTT";

        assert_eq!(actual, expected)
    }

    #[test]
    fn check_ref_amino_ref_matches() {
        let record = PanelRecord {
            gene: "G".to_string(),
            variant: Variant::from_str("C2A").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };
        let refseq = b"AACCTTGTGCAGG".to_vec();
        let padding = 2;

        let actual = record.check_ref(&refseq, padding).unwrap();
        let expected = b"TGT";

        assert_eq!(actual, expected)
    }

    #[test]
    fn check_ref_amino_ref_does_not_match() {
        let record = PanelRecord {
            gene: "G".to_string(),
            variant: Variant::from_str("C2A").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };
        let refseq = b"AACCTTGAGCAGG".to_vec();
        let padding = 2;

        let actual = record.check_ref(&refseq, padding).unwrap_err();
        let expected = PanelError::RefDoesNotMatch(record.name());

        assert_eq!(actual, expected)
    }

    #[test]
    fn check_ref_unknown_amino_ref_does_not_match() {
        let record = PanelRecord {
            gene: "G".to_string(),
            variant: Variant::from_str("Z2A").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };
        let refseq = b"AACCTTGAGCAGG".to_vec();
        let padding = 2;

        let actual = record.check_ref(&refseq, padding).unwrap_err();
        let expected = PanelError::RefDoesNotMatch(record.name());

        assert_eq!(actual, expected)
    }

    #[test]
    fn all_alt_alleles_nucleic_not_x_returns_alt_only() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("A1T").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };

        let actual = record.all_alt_alleles().unwrap();
        let expected = vec![b"T"];
        assert_eq!(actual, expected)
    }

    #[test]
    fn all_alt_alleles_amino_not_x_returns_alt_only() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("A1T").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        let actual = record.all_alt_alleles().unwrap();
        let expected = vec![b"ACT", b"ACC", b"ACA", b"ACG"];
        assert_eq!(actual, expected)
    }

    #[test]
    fn all_alt_alleles_multi_amino_not_x_returns_error() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("A1TT").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        let actual = record.all_alt_alleles().unwrap_err();
        let expected = PanelError::MultiAminoAlleleNotSupported(record.name());
        assert_eq!(actual, expected)
    }

    #[test]
    fn all_alt_alleles_nucleic_x_returns_all_others() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("A1X").unwrap(),
            residue: Residue::Nucleic,
            drugs: Default::default(),
        };

        let actual = record.all_alt_alleles().unwrap();
        let expected = vec![b"C", b"G", b"T"];
        assert_eq!(actual, expected)
    }

    #[test]
    fn all_alt_alleles_amino_x_returns_all_others() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("A1X").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        let mut actual = record.all_alt_alleles().unwrap();
        let mut expected = vec![];
        AMINO_ACIDS
            .iter()
            .filter(|c| *c != &b'A')
            .for_each(|c| expected.extend(amino_to_codons(*c)));
        expected.sort();
        actual.sort();
        assert_eq!(actual, expected)
    }

    #[test]
    fn all_alt_alleles_multi_amino_with_x_returns_error() {
        let record = PanelRecord {
            gene: "".to_string(),
            variant: Variant::from_str("A1RX").unwrap(),
            residue: Residue::Amino,
            drugs: Default::default(),
        };

        let actual = record.all_alt_alleles().unwrap_err();
        let expected = PanelError::MultiAminoAlleleNotSupported(record.name());
        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_is_indel_snp() {
        let v = Variant::from_str("A4T").unwrap();
        assert!(!v.is_indel())
    }

    #[test]
    fn test_variant_is_indel_mnp() {
        let v = Variant::from_str("AA4TA").unwrap();
        assert!(!v.is_indel())
    }

    #[test]
    fn test_variant_is_indel_deletion() {
        let v = Variant::from_str("AA4A").unwrap();
        assert!(v.is_indel())
    }

    #[test]
    fn test_variant_is_indel_insertion() {
        let v = Variant::from_str("AA4ACGT").unwrap();
        assert!(v.is_indel())
    }

    #[test]
    fn test_variant_is_snp_insertion() {
        let v = Variant::from_str("AA4ACGT").unwrap();
        assert!(!v.is_snp())
    }

    #[test]
    fn test_variant_is_snp_snp() {
        let v = Variant::from_str("A4T").unwrap();
        assert!(v.is_snp())
    }

    #[test]
    fn test_variant_is_snp_mnp() {
        let v = Variant::from_str("AA4GT").unwrap();
        assert!(!v.is_snp())
    }

    #[test]
    fn test_variant_simplify_nothing_to_do() {
        let s = "K2*";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = v;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_second_base() {
        let s = "K*2L*";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("K2L").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_first_base() {
        let s = "AR3AK";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("R4K").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_first_two_bases() {
        let s = "CAR3CAK";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("R5K").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_last_two_bases() {
        let s = "CAR3TAR";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("C3T").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_last_two_bases_and_first_two_bases() {
        let s = "QWCAR3QWTAR";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("C5T").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_all_same_does_nothing() {
        let s = "QWCAR3QWCAR";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = v;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_long_tail() {
        let s = "GAGCAG2123CAGCAG";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("G2123C").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_insertion_one_base_ref_does_nothing() {
        let s = "A2AT";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = v;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_insertion_matches_at_start() {
        let s = "AA2AAT";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("A3AT").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_insertion_matches_at_end_and_start() {
        let s = "AA2AATA";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("A3ATA").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_deletion_single_base_alt_does_nothing() {
        let s = "AA2A";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = v;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_deletion_matches_at_end_and_start() {
        let s = "AAT2AT";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("AT3T").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_variant_simplify_deletion_matches_at_end() {
        let s = "AAT2AT";
        let v = Variant::from_str(s).unwrap();

        let actual = v.simplify();
        let expected = Variant::from_str("AT3T").unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_gene_deletion_fmt() {
        let v = Variant::gene_deletion();

        let actual = v.to_string();
        let expected = "gene_absent".to_string();

        assert_eq!(actual, expected)
    }
}
