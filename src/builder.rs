use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::io::{Read, Seek};
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};
use bio::io::{fasta, gff};
use log::{debug, info, warn};
use rust_htslib::bcf;
use structopt::StructOpt;
use thiserror::Error;

use crate::panel::{Panel, PanelError, PanelRecord};
use crate::Runner;
use bio::alphabets::dna::revcomp;

static META: &str = "##";
const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(StructOpt, Debug)]
pub struct Build {
    /// Path to pandora executable. Will try in src/ext or $PATH if not given
    #[structopt(short = "p", long = "pandora", parse(from_os_str), global = true)]
    pub pandora_exec: Option<PathBuf>,
    /// Path to make_prg executable. Will try in src/ext or $PATH if not given
    #[structopt(short = "m", long = "makeprg", parse(from_os_str))]
    pub makeprg_exec: Option<PathBuf>,
    /// Annotation file that will be used to gather information about genes in panel
    #[structopt(short = "a", long = "gff", parse(from_os_str))]
    pub gff_file: PathBuf,
    /// Panel to build index from
    #[structopt(short = "i", long = "panel", parse(from_os_str))]
    pub panel_file: PathBuf,
    /// Reference genome in FASTA format (must be indexed with samtools faidx)
    #[structopt(short = "f", long = "fasta", parse(from_os_str))]
    pub reference_file: PathBuf,
    /// Number of bases of padding to add to start and end of each gene
    #[structopt(short = "P", long = "padding", default_value = "100")]
    pub padding: u32,
    /// Directory to place output
    #[structopt(short = "o", long = "outdir", default_value = ".")]
    pub outdir: PathBuf,
}

/// A collection of custom errors relating to the build component of this package
#[derive(Error, Debug, PartialEq)]
pub enum BuildError {
    /// Contig is missing from the FASTA index
    #[error("Contig {0} is missing from the reference genome index")]
    MissingContig(String),
    /// Couldn't fetch a gene
    #[error("Couldn't fetch contig {0} with start {1} and end {2}")]
    FetchError(String, u64, u64),
    /// GFF record doesn't have a strand
    #[error("No strand information in GFF for {0}")]
    MissingStrand(String),
}

fn load_annotations_for_genes<R>(
    mut reader: gff::Reader<R>,
    genes: &HashSet<String>,
) -> HashMap<String, gff::Record>
where
    R: std::io::Read,
{
    let mut records = reader.records();
    let mut annotations: HashMap<String, gff::Record> = HashMap::new();

    while let Some(Ok(record)) = records.next() {
        if record.feature_type() != "gene" {
            continue;
        }

        if let Some(g) = record.attributes().get("gene") {
            if genes.contains(g) {
                annotations.insert(g.to_string(), record.to_owned());
            }
        }
    }
    annotations
}

fn extract_gene_from_index<R>(
    record: &gff::Record,
    faidx: &mut fasta::IndexedReader<R>,
    padding: u32,
) -> Result<Vec<u8>, BuildError>
where
    R: Read + Seek,
{
    let mut contig_len: u64 = 0;
    for contig in faidx.index.sequences() {
        if contig.name == record.seqname() {
            contig_len = contig.len;
            break;
        }
    }
    if contig_len == 0 {
        return Err(BuildError::MissingContig(record.seqname().to_string()));
    }

    // GFF3 start is 1-based, inclusive. We want to make it 0-based, inclusive
    let start: u64 = max(
        (*record.start() as isize - padding as isize - 1) as isize,
        0,
    ) as u64;
    // GFF3 end is 1-based, inclusive. We want to make it 0-based, exclusive
    let end: u64 = min(record.end() + padding as u64, contig_len);
    if start > end {
        return Err(BuildError::FetchError(
            record.seqname().to_string(),
            start,
            end,
        ));
    }

    if faidx.fetch(record.seqname(), start, end).is_err() {
        return Err(BuildError::FetchError(
            record.seqname().to_string(),
            start,
            end,
        ));
    } else {
        let mut seq: Vec<u8> = Vec::with_capacity((end - start) as usize);
        let is_rev = match record.strand() {
            Some(strand) => strand.strand_symbol() == "-",
            _ => {
                return Err(BuildError::MissingStrand(
                    record.attributes().get("gene").unwrap().to_string(),
                ))
            }
        };
        match faidx.read(&mut seq) {
            Err(_) => Err(BuildError::FetchError(
                record.seqname().to_string(),
                start,
                end,
            )),
            Ok(_) => {
                if is_rev {
                    Ok(revcomp(seq))
                } else {
                    Ok(seq.to_owned())
                }
            }
        }
    }
}

impl Runner for Build {
    fn run(&self) -> Result<()> {
        if !self.outdir.exists() {
            info!("Outdir doesn't exist - creating...");
            std::fs::create_dir(&self.outdir)
                .context(format!("Failed to create {:?}", &self.outdir))?;
        }
        info!("Building panel index...");

        info!("Loading the panel...");
        let mut panel: Panel = HashMap::new();
        let mut genes: HashSet<String> = HashSet::new();
        let mut record_num = 0;

        let mut panel_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(&self.panel_file)?;

        for result in panel_reader.deserialize() {
            record_num += 1;
            let record: PanelRecord = result?;
            genes.insert(record.gene.to_owned());
            let set_of_records = panel
                .entry(record.gene.to_owned())
                .or_insert_with(HashSet::new);
            let seen_before = !set_of_records.insert(record);

            if seen_before {
                warn!(
                    "Duplicate panel record detected in record number {}",
                    record_num
                )
            }
        }
        info!("Loaded {} panel record(s)", panel.len());

        info!("Loading genome annotation for panel genes...");
        debug!("Panel genes: {:?}", genes);
        let gff_reader = gff::Reader::from_file(&self.gff_file, gff::GffType::GFF3)?;
        let annotations = load_annotations_for_genes(gff_reader, &genes);
        if genes.len() != annotations.len() {
            let annotation_genes: HashSet<String> =
                annotations.keys().cloned().collect();
            let diff = genes.difference(&annotation_genes);
            warn!("Couldn't load annotations for genes {:?}", diff);
        }
        info!("Loaded annotations");

        debug!("Creating VCF header...");
        let mut vcf_header = bcf::header::Header::new();
        vcf_header.push_record(format!("{}source=drprgV{}", META, VERSION).as_bytes());
        // add contigs to vcf header
        for (gene, gff_record) in &annotations {
            let length: u64 = (gff_record.end() + 1) - gff_record.start()
                + (&self.padding * 2) as u64;
            vcf_header.push_record(&*vcf_contig_field(gene, length));
        }
        for entry in PanelRecord::vcf_header_entries() {
            vcf_header.push_record(entry);
        }
        debug!("VCF header created");

        let panel_vcf_path = self.outdir.join("panel.bcf");
        let mut vcf_writer = bcf::Writer::from_path(
            &panel_vcf_path,
            &vcf_header,
            false,
            bcf::Format::BCF,
        )?;
        debug!("Loading the reference genome index...");
        let mut faidx = fasta::IndexedReader::from_file(&self.reference_file)?;
        debug!("Loaded the reference genome index");
        let gene_refs_path = self.outdir.join("genes.fa");
        let mut fa_writer = fasta::Writer::to_file(gene_refs_path)?;
        let premsa_dir = self.outdir.join("premsa");
        if !premsa_dir.exists() {
            debug!("Pre-MSA directory doesn't exist - creating...");
            std::fs::create_dir(&premsa_dir)
                .context(format!("Failed to create {:?}", &premsa_dir))?;
        }

        info!("Converting the panel to a VCF...");
        for (gene, gff_record) in &annotations {
            let seq = extract_gene_from_index(&gff_record, &mut faidx, self.padding)?;
            fa_writer.write(gene, Some(&format!("padding={}", self.padding)), &seq)?;
            let premsa_path = premsa_dir.join(format!("{}.premsa.fa", gene));
            let mut premsa_writer = fasta::Writer::to_file(premsa_path)?;
            premsa_writer.write(
                gene,
                Some(&format!("padding={}", self.padding)),
                &seq,
            )?;

            let panel_records_for_gene = &panel[gene];
            for panel_record in panel_records_for_gene {
                let mut vcf_record = vcf_writer.empty_record();
                match panel_record.to_vcf(&mut vcf_record, &seq, self.padding) {
                    Err(PanelError::PosOutOfRange(pos, gene)) => {
                        warn!(
                            "Position {} is out of range for gene {} [Skipping]",
                            pos, gene
                        );
                        continue;
                    }
                    Err(e) => return Err(anyhow!(e)),
                    _ => (),
                }
                // we use unwrap as we have already confirmed above that the strand is + or -
                let strand = gff_record.strand().unwrap().strand_symbol().to_owned();
                vcf_record
                    .push_info_string(b"ST", &[strand.as_bytes()])
                    .context(format!(
                        "Couldn't set INFO field ST for {}",
                        panel_record.name()
                    ))?;
                vcf_writer.write(&vcf_record)?;

                let s: usize = vcf_record.pos() as usize;
                let e: usize = s + vcf_record.inner().rlen as usize;
                for (i, allele) in vcf_record.alleles()[1..].iter().enumerate() {
                    // todo: create mutated seq
                    let mutated_seq =
                        [seq[..s].to_vec(), allele.to_vec(), seq[e..].to_vec()]
                            .concat();
                    // todo: write mutated seq
                    premsa_writer.write(
                        &format!("{}_{}", panel_record.name(), i),
                        None,
                        &mutated_seq,
                    )?;
                }
            }
        }
        info!(
            "Panel successfully converted to a VCF at {:?}",
            panel_vcf_path
        );
        // let makeprg = MakePrg::from_arg(&self.makeprg_exec)?;
        Ok(())
    }
}

fn vcf_contig_field(id: &str, length: u64) -> Vec<u8> {
    format!("{}contig=<ID={},length={}>", META, id, length)
        .as_bytes()
        .to_owned()
}

#[cfg(test)]
mod tests {
    use std::iter::FromIterator;

    use bio::io::fasta::IndexedReader;

    use super::*;

    #[test]
    fn load_annotations_when_no_genes_in_common_returns_empty() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let genes = HashSet::from_iter(vec!["geneX".to_string()]);

        let actual = load_annotations_for_genes(reader, &genes);
        assert!(actual.is_empty())
    }

    #[test]
    fn load_annotations_for_genes_one_gene_in_common() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tgene\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let genes = HashSet::from_iter(vec!["geneX".to_string(), "dnaA".to_string()]);

        let actual = load_annotations_for_genes(reader, &genes);
        let (actual_gene, actual_record) = actual.iter().next().unwrap();
        assert_eq!(actual_gene, "dnaA");
        assert_eq!(actual_record.end(), &1524u64)
    }

    #[test]
    fn load_annotations_for_genes_is_cds_returns_empty() {
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tCDS\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let genes = HashSet::from_iter(vec!["geneX".to_string(), "dnaA".to_string()]);

        let actual = load_annotations_for_genes(reader, &genes);
        assert!(actual.is_empty())
    }

    #[test]
    fn extract_gene_from_index_gene_not_in_index() {
        let padding = 0;
        const GFF: &[u8] = b"NC_000962.3\tRefSeq\tCDS\t1\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap_err();
        let expected = BuildError::MissingContig("NC_000962.3".to_string());
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_interval_out_of_bounds() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t100\t1524\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap_err();
        let expected = BuildError::FetchError("chr1".to_string(), 99, 16);
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_first_base() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t1\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"G";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_too_much_padding_left_wraps_to_start() {
        let padding = 2;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t1\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"GTA";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_too_much_padding_right_wraps_to_end() {
        let padding = 4;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t16\t16\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"ACCCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_padding_start_and_end_exactly_the_same_as_gene() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t16\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"GTAGGCTGAAAACCCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_on_reverse_strand() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t16\t.\t-\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = revcomp(b"GTAGGCTGAAAACCCC");
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_strand() {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t16\t.\t.\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap_err();
        let expected = BuildError::MissingStrand("dnaA".to_string());
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_padding_start_same_as_gene_end_minus_one_from_gene_length(
    ) {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t1\t15\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"GTAGGCTGAAAACCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn extract_gene_from_index_no_padding_start_plus_one_from_gene_start_end_same_as_gene(
    ) {
        let padding = 0;
        const GFF: &[u8] = b"chr1\tRefSeq\tCDS\t2\t16\t.\t+\t.\tID=gene-Rv0001;Dbxref=GeneID:885041;Name=dnaA;experiment=DESCRIPTION:Mutation analysis%2C gene expression[PMID: 10375628];gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=Rv0001\n";
        let mut reader = gff::Reader::new(GFF, gff::GffType::GFF3);
        let record = reader.records().next().unwrap().unwrap();
        const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
        const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
        let mut faidx =
            IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let actual = extract_gene_from_index(&record, &mut faidx, padding).unwrap();
        let expected = b"TAGGCTGAAAACCCC";
        assert_eq!(actual, expected)
    }

    #[test]
    fn test_vcf_contig_field() {
        let id = "chr1";
        let length = 33;

        let actual = vcf_contig_field(id, length);
        let expected = b"##contig=<ID=chr1,length=33>";

        assert_eq!(actual, expected)
    }
}
