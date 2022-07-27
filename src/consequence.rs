use crate::panel::{Residue, Variant};
use crate::report::Evidence;
use bstr::ByteSlice;
use drprg::VcfExt;
use noodles::core::Position;
use noodles::fasta;
use rust_htslib::bcf;
use std::collections::HashMap;

lazy_static! {
    static ref CODON_TABLE: HashMap<&'static [u8], &'static str> = HashMap::from([
        ("TCA".as_bytes(), "S"),
        (b"TCC", "S"),
        (b"TCG", "S"),
        (b"TCT", "S"),
        (b"TTC", "F"),
        (b"TTT", "F"),
        (b"TTA", "L"),
        (b"TTG", "L"),
        (b"TAC", "Y"),
        (b"TAT", "Y"),
        (b"TAA", "*"),
        (b"TAG", "*"),
        (b"TGC", "C"),
        (b"TGT", "C"),
        (b"TGA", "*"),
        (b"TGG", "W"),
        (b"CTA", "L"),
        (b"CTC", "L"),
        (b"CTG", "L"),
        (b"CTT", "L"),
        (b"CCA", "P"),
        (b"CCC", "P"),
        (b"CCG", "P"),
        (b"CCT", "P"),
        (b"CAC", "H"),
        (b"CAT", "H"),
        (b"CAA", "Q"),
        (b"CAG", "Q"),
        (b"CGA", "R"),
        (b"CGC", "R"),
        (b"CGG", "R"),
        (b"CGT", "R"),
        (b"ATA", "I"),
        (b"ATC", "I"),
        (b"ATT", "I"),
        (b"ATG", "M"),
        (b"ACA", "T"),
        (b"ACC", "T"),
        (b"ACG", "T"),
        (b"ACT", "T"),
        (b"AAC", "N"),
        (b"AAT", "N"),
        (b"AAA", "K"),
        (b"AAG", "K"),
        (b"AGC", "S"),
        (b"AGT", "S"),
        (b"AGA", "R"),
        (b"AGG", "R"),
        (b"GTA", "V"),
        (b"GTC", "V"),
        (b"GTG", "V"),
        (b"GTT", "V"),
        (b"GCA", "A"),
        (b"GCC", "A"),
        (b"GCG", "A"),
        (b"GCT", "A"),
        (b"GAC", "D"),
        (b"GAT", "D"),
        (b"GAA", "E"),
        (b"GAG", "E"),
        (b"GGA", "G"),
        (b"GGC", "G"),
        (b"GGG", "G"),
        (b"GGT", "G")
    ]);
}

pub fn consequence_of_variant(
    vcf_record: &bcf::Record,
    padding: i64,
    gene: &fasta::Record,
) -> Result<Evidence, String> {
    if vcf_record.contig() != gene.name() {
        return Err("Contig names don't match".to_string());
    }

    let vcfid = vcf_record.id().to_str_lossy().to_string();
    let ref_allele = vcf_record.alleles()[0];
    let alt_allele = vcf_record.alleles()[vcf_record.called_allele() as usize];
    let is_indel = ref_allele.len() != alt_allele.len();

    // we can unwrap here as VCF doesn't do negative values for POS
    let ref_start = Position::new(vcf_record.pos() as usize + 1).unwrap();
    // unwrap here as it checks for usize overflow, which we are VERY unlikely to ever hit
    let ref_end = ref_start.checked_add(vcf_record.rlen() as usize).unwrap();
    let iv = ref_start..ref_end;
    let gene_seq_at_pos = gene
        .sequence()
        .get(iv)
        .ok_or_else(|| "Could not get gene reference sequence".to_string())?;

    if gene_seq_at_pos != ref_allele {
        return Err(format!(
            "Reference allele {} at position {} doesn't match gene ({}) sequence {}",
            ref_allele.to_str_lossy(),
            vcf_record.pos() + 1,
            gene.name(),
            gene_seq_at_pos.to_str_lossy()
        ));
    }

    // reminder: vcf pos is 0-based. norm_pos is 1-based
    let norm_pos = match vcf_record.pos() {
        p if p < padding => p - padding,
        p => p - (padding - 1),
    };

    let gene_len = gene.sequence().len() as i64 - (padding * 2);
    let var_crosses_gene_end_boundary =
        (norm_pos - 1) + ref_allele.len() as i64 > gene_len;

    if norm_pos.is_negative() || var_crosses_gene_end_boundary || is_indel {
        let variant = Variant {
            reference: ref_allele.to_str_lossy().to_string(),
            pos: norm_pos,
            new: alt_allele.to_str_lossy().to_string(),
        };
        return Ok(Evidence {
            variant: variant.simplify(),
            gene: gene.name().to_string(),
            residue: Residue::Nucleic,
            vcfid,
        });
    }

    // we now know we have a vcf record where the ref and alt are the same length and the variant
    // is located within the coding sequence of the gene
    let gene_start = Position::new((padding + 1) as usize).unwrap();
    let gene_stop = Position::new((gene_len + padding) as usize).unwrap();
    let gene_iv = gene_start..=gene_stop;
    let gene_seq = gene
        .sequence()
        .get(gene_iv)
        .ok_or_else(|| format!("Couldn't extract gene sequence for {}", gene.name()))?;
    let codon_start = (norm_pos - 1) / 3 * 3; // 0-based
    let codon_end = ((norm_pos - 1) + (ref_allele.len() as i64) - 1) / 3 * 3 + 3; // 0-based exclusive
    let codon_iv = codon_start as usize..codon_end as usize;
    let codon_seq = gene_seq
        .get(codon_iv)
        .ok_or("Couldn't extract codon sequence from gene")?;
    let ref_codons = codon_seq.chunks(3);
    let mut mutated_codon_seq = codon_seq.to_owned();
    // norm_pos is 1-based and codon_start is 0-based
    let alt_start_on_codon = (norm_pos - 1) - codon_start;
    mutated_codon_seq.splice(
        alt_start_on_codon as usize..alt_start_on_codon as usize + ref_allele.len(),
        alt_allele.to_vec(),
    );
    let alt_codons = mutated_codon_seq.chunks(3);

    let mut ref_prot = String::new();
    let mut alt_prot = String::new();
    for (r_codon, a_codon) in ref_codons.zip(alt_codons) {
        let r_aa = CODON_TABLE[r_codon];
        ref_prot.push_str(r_aa);
        let a_aa = CODON_TABLE[a_codon];
        alt_prot.push_str(a_aa);
    }

    let codon_num = (norm_pos - 1) / 3 + 1;

    let variant = Variant {
        reference: ref_prot,
        pos: codon_num,
        new: alt_prot,
    };

    Ok(Evidence {
        variant: variant.simplify(),
        gene: gene.name().to_string(),
        residue: Residue::Amino,
        vcfid,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::panel::{Residue, Variant};
    use noodles::fasta::record::{Definition, Sequence};
    use rust_htslib::bcf::record::GenotypeAllele;
    use rust_htslib::bcf::Header;
    use tempfile::NamedTempFile;

    fn generate_pnca() -> fasta::Record {
        let def = Definition::new("pncA", None);
        let seq = Sequence::from(b"GTCATGTTCGCGATCGTCGCGGCGTCATGGACCCTATATCTGTGGCTGCCGCGTCGGTAGGCAAACTGCCCGGGCAGTCGCCCGAACGTATGGTGGACGTATGCGGGCGTTGATCATCGTCGACGTGCAGAACGACTTCTGCGAGGGTGGCTCGCTGGCGGTAACCGGTGGCGCCGCGCTGGCCCGCGCCATCAGCGACTACCTGGCCGAAGCGGCGGACTACCATCACGTCGTGGCAACCAAGGACTTCCACATCGACCCGGGTGACCACTTCTCCGGCACACCGGACTATTCCTCGTCGTGGCCACCGCATTGCGTCAGCGGTACTCCCGGCGCGGACTTCCATCCCAGTCTGGACACGTCGGCAATCGAGGCGGTGTTCTACAAGGGTGCCTACACCGGAGCGTACAGCGGCTTCGAAGGAGTCGACGAGAACGGCACGCCACTGCTGAATTGGCTGCGGCAACGCGGCGTCGATGAGGTCGATGTGGTCGGTATTGCCACCGATCATTGTGTGCGCCAGACGGCCGAGGACGCGGTACGCAATGGCTTGGCCACCAGGGTGCTGGTGGACCTGACAGCGGGTGTGTCGGCCGATACCACCGTCGCCGCGCTGGAGGAGATGCGCACCGCCAGCGTCGAGTTGGTTTGCAGCTCCTGATGGCACCGCCGAACCGGGATGAACTGTTGGCGGCGGTGGAGCGCTCGCCGCAAGCGGCCGCCGCGCACGACCGCGCCGGCTGGGTCGGGTTGTTCACCGG".to_vec());
        fasta::Record::new(def, seq)
    }

    #[test]
    /// Upstream of a gene should just return the evidence as nucleic acid with POS as
    /// a negative (normalised)
    fn test_consequence_of_variant_upstream_of_gene() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "TCA";
        let alt_allele = "TC";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(1);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();

        let expected_variant = Variant {
            reference: "CA".to_string(),
            pos: -98,
            new: "C".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Nucleic,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    /// Upstream of a gene should just return the evidence as nucleic acid with POS as
    /// a negative (normalised). This test checks near the edge of the padding and gene start
    fn test_consequence_of_variant_upstream_of_gene_at_edge() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "TATG";
        let alt_allele = "T";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(99);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();

        let expected_variant = Variant {
            reference: ref_allele.to_string(),
            pos: -1,
            new: alt_allele.to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Nucleic,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    /// Downstream (past stop codon) of a gene should just return the evidence as nucleic acid with POS
    /// normalised.
    fn test_consequence_of_variant_downstream_of_gene() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "CAC";
        let alt_allele = "CA";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(664);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();

        let expected_variant = Variant {
            reference: "AC".to_string(),
            pos: 566,
            new: "A".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Nucleic,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    /// Downstream (past stop codon) of a gene should just return the evidence as nucleic acid with POS
    /// normalised. This checks when the record is the last base after the gene
    fn test_consequence_of_variant_downstream_of_gene_edge() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "TGG";
        let alt_allele = "TG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(661);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();

        let expected_variant = Variant {
            reference: "GG".to_string(),
            pos: 563,
            new: "G".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Nucleic,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_ref_doesnt_match() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "ATG";
        let alt_allele = "GTG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(101);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap_err();
        let expected = format!(
            "Reference allele {} at position {} doesn't match gene ({}) sequence {}",
            ref_allele,
            record.pos() + 1,
            gene.name(),
            "TGC"
        );

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_indel_returned_as_nucleic() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "CGG";
        let alt_allele = "CG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(103);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "GG".to_string(),
            pos: 5,
            new: "G".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Nucleic,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_single_codon_ref_is_whole_codon() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "CGG";
        let alt_allele = "AGG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(103);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "R".to_string(),
            pos: 2,
            new: "R".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_single_base_in_codon() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "G";
        let alt_allele = "T";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(105);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "R".to_string(),
            pos: 2,
            new: "R".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_last_codon() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "TGA";
        let alt_allele = "TTA";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(658);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "*".to_string(),
            pos: 187,
            new: "L".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_last_base() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "A";
        let alt_allele = "T";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(660);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "*".to_string(),
            pos: 187,
            new: "C".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_two_codons_two_changes() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "GCGTTG";
        let alt_allele = "GAGGTG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(106);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "AL".to_string(),
            pos: 3,
            new: "EV".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_two_bases_cross_codon_boundary() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "GA";
        let alt_allele = "CT";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(111);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "LI".to_string(),
            pos: 4,
            new: "FF".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_two_bases_cross_gene_end_boundary() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "AT";
        let alt_allele = "TG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(660);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: ref_allele.to_string(),
            pos: 561,
            new: alt_allele.to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Nucleic,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_five_bases_cross_three_codons() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "GCGGG";
        let alt_allele = "ACAGA";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(102);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "MRA".to_string(),
            pos: 1,
            new: "IQT".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_simplify_multi_codon_synonymous() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "GCGGG";
        let alt_allele = "GCAGA";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(102);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "RA".to_string(),
            pos: 2,
            new: "QT".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_simplify_multi_codon_synonymous_front_and_back() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "GCGGG";
        let alt_allele = "GCAGG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(102);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "R".to_string(),
            pos: 2,
            new: "Q".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_simplify_multi_codon_synonymous_front_two() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "GCGGG";
        let alt_allele = "GCGGA";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(102);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "A".to_string(),
            pos: 3,
            new: "T".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_simplify_two_codon_synonymous_first() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "GCGG";
        let alt_allele = "GTGG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(102);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "R".to_string(),
            pos: 2,
            new: "W".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_simplify_two_codon_synonymous_last() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "GCGG";
        let alt_allele = "ACGG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(102);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "M".to_string(),
            pos: 1,
            new: "I".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_consequence_of_variant_multi_codon_synonymous_all_same_no_simplify() {
        let gene = generate_pnca();
        let vcfid = "id".to_string();
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = Header::new();

        header.push_record(format!("##contig=<ID={}>", gene.name()).as_bytes());
        header.push_sample(b"sample").push_record(
            br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
        );
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();

        let ref_allele = "GCGGG";
        let alt_allele = "GCGTG";

        let alleles: &[&[u8]] = &[ref_allele.as_bytes(), alt_allele.as_bytes()];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record.set_pos(102);
        let rid = record.header().name2rid(gene.name().as_bytes()).unwrap();
        record.set_rid(Some(rid));
        record
            .set_id(vcfid.as_bytes())
            .expect("Failed to set vcf ID");

        let padding = 100;

        let actual = consequence_of_variant(&record, padding, &gene).unwrap();
        let expected_variant = Variant {
            reference: "MRA".to_string(),
            pos: 1,
            new: "MRA".to_string(),
        };
        let expected = Evidence {
            variant: expected_variant,
            gene: gene.name().to_string(),
            residue: Residue::Amino,
            vcfid,
        };

        assert_eq!(actual, expected)
    }
}
