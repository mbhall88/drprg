import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TextIO, Set, Dict, NamedTuple, Optional, List

import click
from cyvcf2 import VCF, Writer
from intervaltree import IntervalTree
from loguru import logger

TRANSLATE = str.maketrans("ATGC", "TACG")


def revcomp(s: str) -> str:
    return s.upper().translate(TRANSLATE)[::-1]


class Genotype(NamedTuple):
    allele1: int
    allele2: int

    def is_null(self) -> bool:
        """Is the genotype null. i.e. ./."""
        return self.allele1 == -1 and self.allele2 == -1

    def is_hom(self) -> bool:
        """Is the genotype homozygous"""
        if self.is_null():
            return False
        if self.allele1 == -1 or self.allele2 == -1:
            return True
        return self.allele1 == self.allele2

    def is_het(self) -> bool:
        """Is the genotype heterozyhous"""
        return not self.is_null() and not self.is_hom()

    def is_hom_ref(self) -> bool:
        """Is genotype homozygous reference?"""
        return self.is_hom() and (self.allele1 == 0 or self.allele2 == 0)

    def is_hom_alt(self) -> bool:
        """Is genotype homozygous alternate?"""
        return self.is_hom() and (self.allele1 > 0 or self.allele2 > 0)

    def alt_index(self) -> Optional[int]:
        """If the genotype is homozygous alternate, returns the 0-based index of the
        alt allele in the alternate allele array.
        """
        if not self.is_hom_alt():
            return None
        return max(self.allele1, self.allele2) - 1

    def allele_index(self) -> Optional[int]:
        """The index of the called allele"""
        if self.is_hom_ref() or self.is_null():
            return 0
        elif self.is_hom_alt():
            return self.alt_index() + 1
        else:
            raise NotImplementedError(f"Het Genotype is unexpected: {self}")

    @staticmethod
    def from_arr(arr: List[int]) -> "Genotype":
        alleles = [a for a in arr if type(a) is int]
        if len(alleles) < 2:
            alleles.append(-1)
        return Genotype(*alleles)


def extract_genes_from_panel(stream: TextIO) -> Set[str]:
    genes = set()
    for line in map(str.rstrip, stream):
        if not line:
            continue
        fields = line.split("\t")
        if gene := fields[0]:
            genes.add(gene)
    return genes


def extract_intervals_for_genes_from_gff(
    genes: Set[str], gff_stream: TextIO, padding: int = 0
) -> IntervalTree:
    intervals = []
    for row in map(str.rstrip, gff_stream):
        if row.startswith("#") or not row:
            continue
        fields = row.split("\t")
        if fields[2].lower() != "gene":
            continue

        attributes = attributes_dict_from_str(fields[8])
        name = attributes.get("gene", attributes.get("Name", None))
        if name is None:
            raise ValueError(f"No gene/Name attribute for ID {attributes['ID']}")
            continue

        if name not in genes:
            continue

        start = (int(fields[3]) - 1) - padding  # GFF start is 1-based inclusive
        end = int(fields[4]) + padding  # GFF end is 1-based inclusive
        strand = fields[6]
        intervals.append((start, end, (name, strand)))

    return IntervalTree.from_tuples(intervals)


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
@click.command()
@click.option(
    "-P",
    "--padding",
    default=100,
    type=int,
    help="Number of bases up and downstream to include for a gene",
    show_default=True,
)
@click.option(
    "-a",
    "--apply-filters",
    is_flag=True,
    help="Whether to skip variants that have a FILTER value that isn't PASS or '.'",
)
@click.option(
    "--only-alt",
    is_flag=True,
    help="Only output variants that have an alternate allele genotype call. The default behaviour is apply all ALTs regardless of whether they are called by any sample",
)
@click.option(
    "--adjust-pos",
    is_flag=True,
    help="Adjust the position of variants that occur in reverse strand genes so they are with respect to the forward strand. This will also reverse complement the alleles",
)
@click.option(
    "-i",
    "--panel",
    help="The panel/catalogue the VCF will be used for. Only genes found in this panel will be extracted.",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "-g",
    "--gff",
    help="The GFF3 annotation for the reference genome",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "-b",
    "--vcf",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The VCF to extract genes from",
)
@click.option(
    "-o",
    "--out-vcf",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The output VCF path",
)
def main(
    padding: int,
    apply_filters: bool,
    only_alt: bool,
    adjust_pos: bool,
    panel: Path,
    gff: Path,
    vcf: Path,
    out_vcf: Path,
):

    logger.info("Extracting gene names from panel...")
    with open(panel) as istream:
        genes = extract_genes_from_panel(istream)

    logger.success(f"Extracted {len(genes)} genes from the panel")

    logger.info("Extracting intervals for genes from GFF...")
    with open(gff) as istream:
        ivtree = extract_intervals_for_genes_from_gff(genes, istream, padding)
    logger.success(f"Intervals extracted for {len(ivtree)} genes")

    logger.info(
        "Extracting those VCF records that fall within the gene intervals and altering "
        "their CHROM and POS accordingly..."
    )
    vcf_reader = VCF(str(vcf))

    for iv in ivtree:
        vcf_reader.add_to_header(f"##contig=<ID={iv.data[0]},length={iv.end-iv.begin}>")

    with TemporaryDirectory() as tmpdirname:
        tmpvcf = str(Path(tmpdirname) / "tmp.vcf")
        vcf_writer = Writer(tmpvcf, tmpl=vcf_reader)

        for record in vcf_reader:
            if apply_filters and record.FILTER is not None:
                continue

            gt = Genotype.from_arr(record.genotypes[0])
            if only_alt and not gt.is_hom_alt():
                continue

            ivs = ivtree[record.start]
            if len(ivs) > 1:
                logger.warning(
                    f"VCF record at POS {record.POS} overlaps with more than 1 gene: {ivs}. "
                    f"Duplicating record - one for each gene..."
                )
            original_record_start = record.start
            original_ref = record.REF
            original_alts = record.ALT
            for iv in ivs:
                chrom, strand = iv.data
                if adjust_pos and strand == "-":
                    norm_pos = (iv.end - original_record_start) - 1
                    ref = revcomp(original_ref)
                    alts = [revcomp(s) for s in original_alts]
                else:
                    norm_pos = original_record_start - iv.begin
                    ref = original_ref
                    alts = original_alts
                record.set_pos(norm_pos)
                record.CHROM = chrom
                record.REF = ref
                record.ALT = alts
                vcf_writer.write_record(record)

        vcf_writer.close()

        logger.info("Sorting VCF...")
        subprocess.run(
            [
                "bcftools",
                "sort",
                "-T",
                tmpdirname,
                "-o",
                str(out_vcf),
                tmpvcf,
            ],
            check=True,
            stderr=sys.stderr,
        )

    vcf_reader.close()

    logger.success("Done!")


if __name__ == "__main__":
    main()
