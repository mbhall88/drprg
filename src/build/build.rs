use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::process::id;

use anyhow::Result;
use bio::io::{fasta, gff};
use log::{debug, info, warn};
use serde::Deserialize;
use structopt::StructOpt;

use drprg::{dependency_dir, MakePrg, MissingDependencies};

use crate::panel::panel::{Panel, PanelRecord};
use crate::Runner;

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

impl Runner for Build {
    fn run(&self) -> Result<()> {
        info!("Building panel index...");

        info!("Loading the panel...");
        let mut panel: Panel = HashSet::new();
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
            let seen_before = !panel.insert(record);

            if seen_before {
                warn!(
                    "Duplicate panel record detected in record number {}",
                    record_num
                )
            }
        }
        info!("Loaded {} panel record(s)", panel.len());

        info!("Loading the reference genome index...");
        let mut faidx = fasta::IndexedReader::from_file(&self.reference_file)?;
        info!("Loaded the reference genome index");

        // open vcf ref file handle
        let gene_refs_path = self.outdir.join("genes.fa");
        let gene_refs_file = std::fs::File::create(gene_refs_path)?;

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
        todo!("Create VCF ref");

        // faidx.fetch_all("NC_000962.3")?;
        // let mut seq: Vec<u8> = vec![];
        // faidx.read(&mut seq);
        // println!("{:?}", seq.into_iter().map(|c| c as char).collect::<String>());

        // let makeprg = MakePrg::from_arg(&self.makeprg_exec)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::iter::FromIterator;

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
}
