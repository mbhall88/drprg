"""This script converts a mykrobe-style panel into the HGVS format accepted by
TBProfiler.
"""
import argparse
import logging
import re
import sys
from dataclasses import dataclass
from difflib import SequenceMatcher
from enum import Enum
from itertools import repeat
from pathlib import Path

BOUNDARY_RGX = re.compile(r"-\d+_\d+")
PROMOTER_DUP_RGX = re.compile(r"-\d+dup[ACGT]")
HEADER = ["Gene", "Mutation", "Drug", "Confers", "Interaction", "Literature"]
STOP = "*"
CODON2AMINO = {
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": STOP,
    "TAG": STOP,
    "TGC": "C",
    "TGT": "C",
    "TGA": STOP,
    "TGG": "W",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
}

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
    STOP: STOP,
}

FRAMESHIFTS: dict[str, str] = snakemake.params.frameshift_genes


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

    @property
    def ref(self) -> str:
        return self.mutation.ref

    @property
    def alt(self) -> str:
        return self.mutation.alt

    @property
    def pos(self) -> int:
        return self.mutation.pos

    def convert(self, biotype: BioType) -> "HgvsVariant":
        return HgvsVariant.from_mykrobe(self, biotype)

    def is_snp(self) -> bool:
        return (
            len(self.ref) == len(self.alt)
            and len(self.ref) == 1
            and self.residue is Residue.Nucleic
        )

    def is_mnp(self) -> bool:
        return (
            len(self.ref) == len(self.alt)
            and len(self.ref) > 1
            and self.residue is Residue.Nucleic
        )

    def is_indel(self) -> bool:
        return len(self.ref) != len(self.alt) and self.residue is Residue.Nucleic

    def is_frameshift(self) -> bool:
        return self.residue is Residue.Nucleic and (
            abs(len(self.ref) - len(self.alt)) % 3 != 0
        )

    def is_ambiguous(self) -> bool:
        return "X" in self.alt

    def disambiguate(self) -> list["MykrobeVariant"]:
        if not self.is_ambiguous():
            return [self]

        if self.residue is Residue.Protein:
            alts = [
                aa for aa in PROTEIN_LETTERS_1TO3.keys() if aa not in (self.ref, STOP)
            ]
        else:
            alts = list("ACGT".replace(self.ref, ""))

        variants = []
        for alt in alts:
            mut = MykrobeMutation(self.ref, self.pos, alt)
            variants.append(MykrobeVariant(self.gene, mut, self.residue))

        return variants


class IndelType(Enum):
    Deletion = "del"
    Duplication = "dup"
    Insertion = "ins"

    @staticmethod
    def from_mutation(mutation: MykrobeMutation) -> "IndelType":
        if len(mutation.ref) > len(mutation.alt):
            return IndelType.Deletion
        elif len(mutation.ref) < len(mutation.alt):
            if (
                set(mutation.ref) == set(mutation.alt)
                and len(mutation.alt) == 2
                and len(mutation.ref) == 1
            ):
                return IndelType.Duplication
            else:
                return IndelType.Insertion
        else:
            raise ValueError(f"{mutation} is not an indel")


@dataclass(frozen=True, eq=True)
class HgvsVariant:
    gene: str
    mutation: str

    def spans_gene_boundary(self) -> bool:
        m = BOUNDARY_RGX.search(self.mutation)
        return m is not None

    @staticmethod
    def from_mykrobe(variant: MykrobeVariant, biotype: BioType) -> "HgvsVariant":
        if variant.residue is Residue.Protein:
            prefix = "p."
            ref = PROTEIN_LETTERS_1TO3[variant.ref]
            alt = PROTEIN_LETTERS_1TO3[variant.alt]
            pos = variant.mutation.pos
            if alt != "*" and len(ref) != len(alt):
                raise NotImplementedError(
                    f"Cannot handle non-substitution protein mutations: {variant.mutation}"
                )
            mut = f"{ref}{pos}{alt}"
        else:
            prefix = "n." if biotype is not BioType.Coding else "c."
            ref = variant.ref
            alt = variant.alt
            pos = variant.pos
            if variant.is_snp():
                mut = f"{pos}{ref}>{alt}"
            elif variant.is_mnp():

                if (pos - 1) % 3 != 0 or len(ref) != 3:
                    raise NotImplementedError(
                        f"Don't know how to deal with MNPs that aren't at the start of a codon: {variant}"
                    )
                codon = pos2codon(pos)
                ref_aa = PROTEIN_LETTERS_1TO3[CODON2AMINO[ref]]
                alt_aa = PROTEIN_LETTERS_1TO3[CODON2AMINO[alt]]
                mut = f"{ref_aa}{codon}{alt_aa}"
                prefix = "p."
            elif variant.is_indel():
                indel_type = IndelType.from_mutation(variant.mutation)
                match indel_type:
                    case IndelType.Insertion:
                        sm = SequenceMatcher(a=ref, b=alt)
                        insertion = next(
                            (tup for tup in sm.get_opcodes() if tup[0] == "insert"),
                            None,
                        )

                        if insertion is None:
                            raise ValueError(f"Expected an insertion but got {variant}")

                        ins_start = pos + (insertion[1] - 1)
                        ins_end = ins_start + 2 if ins_start == -1 else ins_start + 1

                        ins_seq = alt[insertion[3] : insertion[4]]
                        mut = (
                            f"{ins_start}_{ins_end}{IndelType.Insertion.value}{ins_seq}"
                        )
                    case IndelType.Duplication:
                        mut = f"{pos}{IndelType.Duplication.value}{ref}"
                    case IndelType.Deletion:
                        sm = SequenceMatcher(a=ref, b=alt)
                        # we make an assumption here that there is only one deletion block between ref and alt
                        deletion = next(
                            (tup for tup in sm.get_opcodes() if tup[0] == "delete"),
                            None,
                        )

                        if deletion is None:
                            raise ValueError(f"Expected a deletion but got {variant}")

                        del_size = len(ref) - len(alt)
                        del_start = pos + deletion[1]
                        del_end = del_start + (
                            del_size - 1
                        )  # end is inclusive for hgvs
                        # if deletion crosses gene boundary, we need to change the offset
                        if 0 in list(range(del_start, del_end)):
                            if del_start >= 0:
                                del_start += 1
                            del_end += 1

                        del_range = (
                            str(del_start)
                            if del_end == del_start
                            else f"{del_start}_{del_end}"
                        )
                        mut = f"{del_range}{IndelType.Deletion.value}"
                    case _:
                        raise NotImplementedError(variant)

            else:
                raise NotImplementedError(variant)

        return HgvsVariant(gene=variant.gene, mutation=f"{prefix}{mut}")


def pos2codon(pos: int) -> int:
    """Convert a (one-based) genomic position to a (one-based) codon position"""
    pos -= 1
    codon = pos // 3
    return codon + 1


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
    skipped_boundary = 0
    skipped_prom_dup = 0
    converted = dict()
    with open(args.output, "w") as out_fp, open(args.panel) as in_fp:
        print(",".join(HEADER), file=out_fp)

        # add framshift genes
        for gene, drug in FRAMESHIFTS:
            row = [gene, "frameshift", drug, "resistance", "", ""]
            print(",".join(row), file=out_fp)
            counter += 1

        for row in map(str.strip, in_fp):
            if counter % 1000 == 0:
                logging.debug(f"Converted {counter} variants...")

            fields = row.split("\t")
            drug = fields[-1]
            mykrobe_var = MykrobeVariant(
                gene=fields[0],
                mutation=MykrobeMutation.from_str(fields[1]),
                residue=Residue(fields[2]),
            )

            # skip frameshifts as we add them at the gene level
            if mykrobe_var.is_frameshift() and mykrobe_var.gene in FRAMESHIFTS:
                continue

            biotype = biotypes.get(mykrobe_var.gene)
            if biotype is None:
                raise KeyError(f"Could not find a gene biotype for {mykrobe_var.gene}")

            variants = mykrobe_var.disambiguate()

            for var in variants:
                try:
                    hgvs_var = converted.get(var, var.convert(biotype))
                except ValueError as err:
                    logging.warning(f"Failed to convert {var}..skipping..\n{err}")
                    continue

                if hgvs_var.spans_gene_boundary():
                    logging.debug(
                        f"Skipping {var} ({hgvs_var}) as it spans a gene boundary and "
                        f"tb-profiler can't deal with that"
                    )
                    skipped_boundary += 1
                    continue

                is_promoter_duplication = PROMOTER_DUP_RGX.search(hgvs_var.mutation)
                if is_promoter_duplication:
                    logging.debug(
                        f"Skipping {var} ({hgvs_var}) tb-profiler can't deal with "
                        f"duplications in promoters"
                    )
                    skipped_prom_dup += 1
                    continue

                converted[var] = hgvs_var
                row = [hgvs_var.gene, hgvs_var.mutation, drug, "resistance", "", ""]
                print(",".join(row), file=out_fp)
                counter += 1

    logging.info(f"Finished converting {counter} variants")
    logging.info(f"Skipped {skipped_boundary} variants as they spanned a gene boundary")
    logging.info(
        f"Skipped {skipped_prom_dup} variants as they are promoter duplications"
    )


def test():
    test_name = "Amino acid substitution"
    mut = MykrobeMutation.from_str("S450L")
    gene = "rpoB"
    residue = Residue.Protein
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "p.Ser450Leu")

    assert actual == expected, (test_name, actual, expected)

    test_name = "Nucelic acid SNP in coding region"
    mut = MykrobeMutation.from_str("T196G")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.196T>G")

    assert actual == expected, (test_name, actual, expected)

    test_name = "Nucelic acid MNP in coding region"
    mut = MykrobeMutation.from_str("TCG196TAG")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "p.Ser66*")

    assert actual == expected, (test_name, actual, expected)

    test_name = "Nucelic acid SNP in non-coding gene"
    mut = MykrobeMutation.from_str("A1325C")
    gene = "rrs"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.RibosomalRNA

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "n.1325A>C")

    assert actual == expected, (test_name, actual, expected)

    # test deletions
    test_name = "1bp deletion - first base is deleted"
    mut = MykrobeMutation.from_str("AC1325C")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1325del")

    assert actual == expected, (test_name, actual, expected)

    test_name = "1bp deletion - second base is deleted"
    mut = MykrobeMutation.from_str("AC1325A")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1326del")

    assert actual == expected, (test_name, actual, expected)

    test_name = "2bp deletion - first two bases deleted"
    mut = MykrobeMutation.from_str("ACT1325T")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1325_1326del")

    assert actual == expected, (test_name, actual, expected)

    test_name = "1bp deletion - three base ref, first base deleted"
    mut = MykrobeMutation.from_str("TGG884GG")
    gene = "ethA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.884del")

    assert actual == expected, (test_name, actual, expected)

    test_name = "1bp deletion - three base ref, last base deleted"
    mut = MykrobeMutation.from_str("GTA714GT")
    gene = "tylA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.716del")

    assert actual == expected, (test_name, actual, expected)

    test_name = "1bp deletion - three base ref, deletion in the middle"
    mut = MykrobeMutation.from_str("GTA825GA")
    gene = "ethA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.826del")

    assert actual == expected, (test_name, actual, expected)

    test_name = "2bp deletion - four base ref, deletion in the middle"
    mut = MykrobeMutation.from_str("GTAA825GA")
    gene = "ethA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.826_827del")

    assert actual == expected, (test_name, actual, expected)

    test_name = "2bp deletion - second two bases deleted"
    mut = MykrobeMutation.from_str("ACT1325A")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1326_1327del")

    assert actual == expected, (test_name, actual, expected)

    # test insertions (and duplications)
    test_name = "1bp insertion"
    mut = MykrobeMutation.from_str("A1325AC")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1325_1326insC")

    assert actual == expected, (test_name, actual, expected)

    test_name = "1bp insertion - duplication"
    mut = MykrobeMutation.from_str("A1325AA")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1325dupA")

    assert actual == expected, (test_name, actual, expected)

    test_name = "2bp insertion"
    mut = MykrobeMutation.from_str("A1325ACT")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1325_1326insCT")

    assert actual == expected, (test_name, actual, expected)

    test_name = "2bp insertion - duplication"
    mut = MykrobeMutation.from_str("A1325AAA")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1325_1326insAA")

    assert actual == expected, (test_name, actual, expected)

    test_name = "1bp insertion - reference is two bases"
    mut = MykrobeMutation.from_str("GC1112GGC")
    gene = "ethA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1111_1112insG")

    assert actual == expected, (test_name, actual, expected)

    # test promoter
    test_name = "2bp insertion - promoter"
    mut = MykrobeMutation.from_str("A-1325AAA")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.-1325_-1324insAA")

    assert actual == expected, (test_name, actual, expected)

    test_name = "2bp insertion - promoter at gene boundary"
    mut = MykrobeMutation.from_str("A-1ATG")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.-1_1insTG")

    assert actual == expected, (test_name, actual, expected)

    test_name = "2bp deletion - promoter at gene boundary"
    mut = MykrobeMutation.from_str("AGT-1A")
    gene = "pncA"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.1_2del")

    assert actual == expected, (test_name, actual, expected)

    test_name = "Substitution in promoter"
    mut = MykrobeMutation.from_str("C-15T")
    gene = "fabG1"
    residue = Residue.Nucleic
    var = MykrobeVariant(gene, mut, residue)
    biotype = BioType.Coding

    actual = var.convert(biotype)
    expected = HgvsVariant(gene, "c.-15C>T")

    assert actual == expected, (test_name, actual, expected)

    print("All tests pass")


if __name__ == "__main__":
    main()
