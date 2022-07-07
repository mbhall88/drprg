import sys

sys.stderr = open(snakemake.log[0], "w")
from typing import TextIO, Set, Dict
from pysam import FastaFile
from loguru import logger
import subprocess


def extract_genes_from_panel(stream: TextIO) -> Set[str]:
    genes = set()
    for line in map(str.rstrip, stream):
        if not line:
            continue
        fields = line.split("\t")
        if gene := fields[0]:
            genes.add(gene)
    return genes


def attributes_dict_from_str(s: str) -> Dict[str, str]:
    d = dict()
    for pairs in s.split(";"):
        k, v = pairs.split("=")
        if k in d:
            raise KeyError(f"Attribute key {k} appears twice")
        d[k] = v
    return d


##########################################################
# MAIN
##########################################################
def main():
    padding: int = snakemake.params.padding

    logger.info("Extracting gene names from panel...")
    with open(snakemake.input.panel) as istream:
        genes = extract_genes_from_panel(istream)

    logger.success(f"Extracted {len(genes)} genes from the panel")

    logger.info("Writing reference sequences for each gene...")
    faidx = FastaFile(snakemake.input.genome)
    with open(snakemake.input.annotation) as istream, open(
        snakemake.output.fasta, "w"
    ) as ostream:
        for row in map(str.rstrip, istream):
            if row.startswith("#") or not row:
                continue
            fields = row.split("\t")
            if fields[2].lower() != "gene":
                continue

            attributes = attributes_dict_from_str(fields[8])
            name = attributes.get("gene", attributes.get("Name", None))
            if name is None:
                logger.warning(f"No gene/Name attribute for ID {attributes['ID']}")
                continue

            if name not in genes:
                continue

            strand = fields[6]
            chrom = fields[0]
            start = (int(fields[3]) - 1) - padding  # GFF start is 1-based inclusive
            end = int(fields[4]) + padding  # GFF end is 1-based inclusive
            seq = faidx.fetch(reference=chrom, start=start, end=end)
            header = f">{name} strand={strand} chrom={chrom} start={start} end={end}"
            print(f"{header}\n{seq}", file=ostream)
            logger.debug("Wrote {} to file", name)

    logger.info("Indexing fasta...")
    subprocess.run(
        [
            "samtools",
            "faidx",
            snakemake.output.fasta,
        ],
        check=True,
        stderr=sys.stderr,
    )

    logger.success("Done!")


main()
