import re
import sys
from dataclasses import dataclass
from enum import Enum
from itertools import repeat
from typing import Dict, Tuple, TextIO, List

sys.stderr = open(snakemake.log[0], "w")

TRANSLATE = str.maketrans("ATGC", "TACG")
STOP = "*"
Contig = str
Seq = str
Index = Dict[Contig, Seq]


def eprint(msg: str):
    print(msg, file=sys.stderr)


class DuplicateContigsError(Exception):
    pass


codon2amino = {
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


def revcomp(s: str) -> str:
    return complement(s)[::-1]


def complement(s: str) -> str:
    return s.upper().translate(TRANSLATE)


class Strand(Enum):
    Forward = "+"
    Reverse = "-"
    NotRelevant = "."
    Unknown = "?"

    def __str__(self) -> str:
        return str(self.value)


def translate(seq: str, stop_last=True) -> str:
    if len(seq) % 3 != 0:
        raise ValueError("Sequence length must be a multiple of 3")

    prot = ""
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        prot += codon2amino[codon]

    if stop_last and not prot.endswith(STOP):
        raise ValueError("Sequence did not end in a stop codon")

    return prot


@dataclass
class GffFeature:
    seqid: Contig
    source: str
    method: str  # correct term is type, but that is a python reserved variable name
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    score: float
    strand: Strand
    phase: int
    attributes: Dict[str, str]

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
            strand=Strand(fields[6]),
            phase=phase,
            attributes=attributes,
        )

    def slice(self, zero_based: bool = True) -> Tuple[int, int]:
        """Get a tuple for slicing a python object.
        The reason this method is required is that GFF uses 1-based INCLUSIVE
        coordinates. Meaning the end position is also included in the slice.
        """
        if zero_based:
            return self.start - 1, self.end
        return self.start, self.end + 1

    def _extract_sequence(
        self, index: Index, start_offset: int = 0, end_offset: int = 0
    ) -> str:
        refseq = index.get(self.seqid)
        if refseq is None:
            raise IndexError(f"Contig {self.seqid} does not exist in reference")
        s, e = self.slice(zero_based=True)
        s -= start_offset
        e += end_offset
        return refseq[s:e]

    def nucleotide_sequence(
        self, index: Index, start_offset: int = 0, end_offset: int = 0
    ) -> str:
        nuc_seq = self._extract_sequence(
            index, start_offset=start_offset, end_offset=end_offset
        )
        if self.strand is Strand.Reverse:
            nuc_seq = revcomp(nuc_seq)

        return nuc_seq

    def protein_sequence(self, index: Index) -> str:
        nuc_seq = self.nucleotide_sequence(index)
        return translate(nuc_seq)

    @property
    def name(self) -> str:
        for key in ["Name", "gene", "ID"]:
            if name := self.attributes.get(key, ""):
                break

        return name


def is_variant_valid(
    gene: str, variant: str, alphabet: str, index: Index, feature: GffFeature
):
    ref, original_pos, alt = split_var_name(variant)
    if original_pos >= 1:
        pos = original_pos - 1
    else:
        pos = original_pos

    if alphabet == "PROT":
        refseq = feature.protein_sequence(index)
        expected_ref = refseq[pos]
        if ref != expected_ref:
            eprint(
                f"[WARN]: Variant amino acid {ref} does not match reference {expected_ref} at position {original_pos} in {gene}"
            )
            return False
        else:
            return True

    offset = 0 if pos >= 0 else pos - len(ref)
    refseq = feature.nucleotide_sequence(index, start_offset=offset, end_offset=offset)

    if pos < 0:
        expected_ref = refseq[: len(ref)]
    else:
        expected_ref = refseq[pos : pos + len(ref)]

    if ref != expected_ref:
        eprint(
            f"[WARN]: Variant nucleic acid {ref} does not match reference {expected_ref} at position {original_pos} in {gene}"
        )
        return False
    else:
        return True


def index_fasta(stream: TextIO) -> Index:
    fasta_index: Index = dict()
    sequence: List[Seq] = []
    name: Contig = ""
    for line in map(str.rstrip, stream):
        if not line:
            continue
        if line.startswith(">"):
            if sequence and name:
                fasta_index[name] = "".join(sequence)
                sequence = []
            name = line.split()[0][1:]
            if name in fasta_index:
                raise DuplicateContigsError(
                    f"Contig {name} occurs multiple times in the fasta file."
                )
            continue
        else:
            sequence.append(line)
    if name and sequence:
        fasta_index[name] = "".join(sequence)

    return fasta_index


def split_var_name(name: str) -> Tuple[str, int, str]:
    items = re.match(r"([A-Z]+)([-\d]+)([A-Z/*]+)", name, re.I).groups()
    return items[0], int(items[1]), items[2]


def main():
    exclude = set(snakemake.params.exclude)
    out_fp = open(snakemake.output.panel, "w")

    panel_genes: set[str] = set()
    with open(snakemake.input.panel) as fp:
        for line in fp:
            gene = line.split("\t")[0]
            panel_genes.add(gene)
            out_fp.write(line)

    with open(snakemake.input.reference) as fp:
        index = index_fasta(fp)

    features = dict()
    with open(snakemake.input.features) as fp:
        for line in map(str.rstrip, fp):
            if not line or line.startswith("#"):
                continue

            feature = GffFeature.from_str(line)
            if feature.method != "gene":
                continue

            if not feature.name or feature.name not in panel_genes:
                continue

            features[feature.name] = feature

    eprint(f"Loaded {len(panel_genes)} genes from panel")
    c = 0
    with open(snakemake.input.known) as fp:
        _ = next(fp)  # skip header
        for line in map(str.rstrip, fp):
            gene, var, alpha, drug, grading = line.split("\t")
            grading = int(grading)
            mut = f"{gene}_{var}"

            if gene not in panel_genes or mut in exclude:
                continue

            if grading in snakemake.params.keep_grades:
                # The position of DNA mutations (inside genes) on the rev strand point to the *end* of the
                # reference allele. See the data cleaning notebook for an explanation
                is_promotor_mut = "-" in var
                if (
                    features[gene].strand is Strand.Reverse
                    and alpha == "DNA"
                    and not is_promotor_mut
                ):
                    ref, pos, alt = split_var_name(var)
                    ref = revcomp(ref)
                    alt = revcomp(alt)
                    pos -= len(ref) - 1
                    var = f"{ref}{pos}{alt}"

                if is_variant_valid(gene, var, alpha, index, features[gene]):
                    print(
                        "\t".join([gene, var, alpha, snakemake.params.no_drug]),
                        file=out_fp,
                    )
                    c += 1
                else:
                    continue

    eprint(f"Added {c} non-resistant mutations to panel")

    out_fp.close()


main()
