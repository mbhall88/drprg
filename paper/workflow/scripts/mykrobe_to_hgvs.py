"""This script converts a mykrobe-style panel into the HGVS format accepted by
TBProfiler.
"""
import re
import sys
import argparse
from dataclasses import dataclass
from enum import Enum
from itertools import repeat
from pathlib import Path
import logging

PROTEIN_LETTERS_1TO3 = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
    "*": "*",
}


class BioType(Enum):
    Other = "other"
    NonCodingRNA = "ncRNA"
    RibosomalRNA = "rRNA"
    Coding = "protein_coding"
    MiscRNA = "misc_RNA"
    TransferRNA = "tRNA"


@dataclass
class GffFeature:
    seqid: str
    source: str
    method: str  # correct term is type, but that is a python reserved variable name
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    score: float
    strand: str
    phase: int
    attributes: dict[str, str]

    @property
    def name(self) -> str:
        for key in ["Name", "gene", "ID"]:
            if name := self.attributes.get(key, ""):
                break

        return name

    @staticmethod
    def from_str(s: str) -> "GffFeature":
        fields = s.split("\t")
        score = 0 if fields[5] == "." else float(fields[5])
        phase = -1 if fields[7] == "." else int(fields[7])
        attr_fields = fields[-1].split(";")
        attributes = {k: v for k, v in map(str.split, attr_fields, repeat("="))}
        return GffFeature(
            seqid=fields[0],
            source=fields[1],
            method=fields[2],
            start=int(fields[3]),
            end=int(fields[4]),
            score=score,
            strand=fields[6],
            phase=phase,
            attributes=attributes,
        )


def load_biotypes(gff: Path) -> dict[str, BioType]:
    biotypes = dict()
    with open(gff) as fp:
        for line in map(str.strip, fp):
            if not line or line.startswith("#"):
                continue
            ftr = GffFeature.from_str(line)
            if ftr.method != "gene":
                continue

            biotype = ftr.attributes.get("gene_biotype")
            if biotype is None:
                logging.warning(f"No gene biotype found for {ftr.name}")
                continue
            assert ftr.name not in biotypes, ftr
            biotypes[ftr.name] = BioType(biotype)

    return biotypes


class Residue(Enum):
    Protein = "PROT"
    Nucleic = "DNA"


@dataclass(frozen=True)
class MykrobeMutation:
    ref: str
    pos: int
    alt: str

    def __str__(self):
        return f"{self.ref}{self.pos}{self.alt}"

    @staticmethod
    def from_str(s: str) -> "MykrobeMutation":
        items = re.match(r"([A-Z]+)([-\d]+)([A-Z/*]+)", s, re.I).groups()
        return MykrobeMutation(ref=items[0], alt=items[2], pos=int(items[1]))


@dataclass(frozen=True, eq=True)
class MykrobeVariant:
    gene: str
    mutation: MykrobeMutation
    residue: Residue

    def convert(self, biotype: BioType) -> "HgvsVariant":
        if biotype is not BioType.Coding:
            raise NotImplementedError(biotype)
        if self.residue is Residue.Protein:
            prefix = "p."
            ref = PROTEIN_LETTERS_1TO3[self.mutation.ref]
            alt = PROTEIN_LETTERS_1TO3[self.mutation.alt]
            pos = self.mutation.pos
            if alt != "*" and len(ref) != len(alt):
                raise NotImplementedError(
                    f"Cannot handle non-substitution protein mutations: {self.mutation}"
                )
            mut = f"{prefix}{ref}{pos}{alt}"
        else:
            raise NotImplementedError(self.residue)

        return HgvsVariant(gene=self.gene, mutation=mut)


@dataclass(frozen=True, eq=True)
class HgvsVariant:
    gene: str
    mutation: str


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--panel",
        required=True,
        type=Path,
        help=(
            "The (tab-delimited) mykrobe-style panel to convert. This must contain four "
            "columns: gene, mutation, alphabet, drug. "
            "e.g., pncA\tT142R\tPROT\tPyrazinamide"
        ),
    )
    parser.add_argument(
        "-g",
        "--gff",
        type=Path,
        required=True,
        help="GFF file for the panel's reference genome",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="File to write hgvs panel to [default: stdout]",
        type=Path,
    )
    parser.add_argument("--test", help=argparse.SUPPRESS, action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if args.test:
        test()
        sys.exit(0)

    log_lvl = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        format="[%(levelname)s - %(asctime)s]: %(message)s", level=log_lvl
    )

    logging.info("Loading gene biotypes...")
    biotypes = load_biotypes(args.gff)
    logging.info(f"Loaded {len(biotypes)} gene biotypes")

    if args.output is None:
        args.output = 1

    logging.info("Converting the panel...")
    counter = 0
    converted = dict()
    with open(args.output, "w") as out_fp, open(args.panel) as in_fp:
        for row in map(str.strip, in_fp):
            counter += 1
            if counter % 1000 == 0:
                logging.debug(f"Converted {counter} variants...")

            fields = row.split("\t")
            mykrobe_var = MykrobeVariant(
                gene=fields[0],
                mutation=MykrobeMutation.from_str(fields[1]),
                residue=Residue(fields[2]),
            )
            biotype = biotypes.get(mykrobe_var.gene)
            if biotype is None:
                raise KeyError(f"Could not find a gene biotype for {mykrobe_var.gene}")

            hgvs_var = converted.get(mykrobe_var, mykrobe_var.convert(biotype))
            converted[mykrobe_var] = hgvs_var

            drug = fields[-1]

    logging.info(f"Finished converting {counter} variants")


def test():
    test_name = "Amino acid substitution"
    mut = MykrobeMutation.from_str("S450L")
    gene = "rpoB"
    residue = Residue.Protein
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "p.Ser450Leu")

    assert actual == expected, test_name

    test_name = "Nucelic acid SNP in coding region"
    mut = MykrobeMutation.from_str("T196G")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.196T>G")

    assert actual == expected, test_name

    test_name = "Nucelic acid MNP in coding region"
    mut = MykrobeMutation.from_str("TCG196TAG")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.196_198delinsTAG")

    assert actual == expected, test_name

if __name__ == "__main__":
    main()
