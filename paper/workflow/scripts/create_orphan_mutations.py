import sys

sys.stderr = open(snakemake.log[0], "w")
import re
import shlex
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from itertools import repeat
from pathlib import Path
from tempfile import TemporaryDirectory
import difflib


REVERSE = "-"
MISSING = "."
TRANSLATE = str.maketrans("ATGC", "TACG")
STOP = "*"
CHROM = "NC_000962.3"
CODONTAB = {
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

VCF_METALINES = """##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=NC_000962.3,length=4411532>"""

VCF_HEADER = "\t".join(
    ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
)

AMINOTAB = defaultdict(list)
for _codon, _aa in CODONTAB.items():
    AMINOTAB[_aa].append(_codon)


def eprint(msg: str):
    print(msg, file=sys.stderr)


def revcomp(s: str) -> str:
    return complement(s)[::-1]


def complement(s: str) -> str:
    return s.upper().translate(TRANSLATE)


def translate(seq: str, stop_last=True) -> str:
    if len(seq) % 3 != 0:
        raise ValueError("Sequence length must be a multiple of 3")

    prot = ""
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        prot += CODONTAB[codon]

    if stop_last and not prot.endswith(STOP):
        raise ValueError("Sequence did not end in a stop codon")

    return prot


def get_closest_codon(from_codon: str, to_aa: str) -> str:
    return difflib.get_close_matches(from_codon, AMINOTAB[to_aa], n=1, cutoff=0)[0]


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

    def slice(self, zero_based: bool = True) -> tuple[int, int]:
        """Get a tuple for slicing a python object.
        The reason this method is required is that GFF uses 1-based INCLUSIVE
        coordinates. Meaning the end position is also included in the slice.
        """
        if zero_based:
            return self.start - 1, self.end
        return self.start, self.end + 1

    def _extract_sequence(
        self, index: dict[str, str], start_offset: int = 0, end_offset: int = 0
    ) -> str:
        refseq = index.get(self.seqid)
        if refseq is None:
            raise IndexError(f"Contig {self.seqid} does not exist in reference")
        s, e = self.slice(zero_based=True)
        s -= start_offset
        e += end_offset
        return refseq[s:e]

    def nucleotide_sequence(
        self, index: dict[str, str], start_offset: int = 0, end_offset: int = 0
    ) -> str:
        nuc_seq = self._extract_sequence(
            index, start_offset=start_offset, end_offset=end_offset
        )
        if self.strand == REVERSE:
            nuc_seq = revcomp(nuc_seq)

        return nuc_seq

    def protein_sequence(self, index: dict[str, str]) -> str:
        nuc_seq = self.nucleotide_sequence(index)
        return translate(nuc_seq)

    @property
    def name(self) -> str:
        for key in ["Name", "gene", "ID"]:
            if name := self.attributes.get(key, ""):
                break

        return name


def split_var_name(name: str) -> tuple[str, int, str]:
    if "ins" in name or "del" in name:
        items = name.split("_")
        return items[3], int(items[0]), items[-1]

    items = re.match(r"([A-Z]+)([-\d]+)([A-Z/*]+)", name, re.I).groups()
    return items[0], int(items[1]), items[2]


def main():
    reference = ""
    with open(snakemake.input.reference) as fp:
        for line in map(str.rstrip, fp):
            if line and line[0] == ">":
                continue
            reference += line

    features = dict()
    with open(snakemake.input.annotation) as fp:
        for line in map(str.rstrip, fp):
            if not line or line.startswith("#"):
                continue

            feature = GffFeature.from_str(line)
            if feature.method != "gene":
                continue

            features[feature.name] = feature

    tmpdir = TemporaryDirectory()
    tmpdirname = Path(tmpdir.name)

    vcfs_to_merge = []

    with open(snakemake.input.mutations) as fp:
        for line in map(str.rstrip, fp):
            eprint(f"Processing {line}...")
            gene, mut = line.split("_", maxsplit=1)
            mut = mut.replace("!", STOP)
            ref, pos, alt = split_var_name(mut)
            is_dna = ref.islower()
            ref = ref.upper()
            alt = alt.upper()
            ftr = features[gene]
            is_indel = len(ref) != len(alt)

            if is_dna:
                mut_gene_pos = pos + 1 if pos < 1 else pos
            else:
                mut_gene_pos = pos * 3 - 2

            is_reverse_strand = ftr.strand == REVERSE
            l = len(ref) if is_dna else 3
            if is_reverse_strand:
                _start = ftr.end - mut_gene_pos
                if is_indel:
                    _end = _start + l
                    start = _start
                    end = _end
                else:
                    _end = _start - l
                    start = _end + 1
                    end = _start + 1
            else:
                start = ftr.start + (mut_gene_pos - 1)
                start -= 1
                end = start + l

            refseq = reference[start:end]
            if is_reverse_strand and not is_indel:
                refseq = revcomp(refseq)
            if is_dna:
                if ref != refseq:
                    raise ValueError(
                        f"{line}: ref allele {ref} does not match reference {refseq}"
                    )
            else:
                aa = CODONTAB[refseq]
                if aa != ref:
                    raise ValueError(
                        f"{line}: ref allele {ref} does not match reference amino acid {aa}"
                    )

            vcf_ref = reference[start:end]
            pos = start + 1
            if is_indel:
                assert vcf_ref == ref
                vcf_alt = alt
            elif is_dna:
                vcf_alt = revcomp(alt) if is_reverse_strand else alt
            else:
                vcf_alt = get_closest_codon(refseq, alt)
                if is_reverse_strand:
                    vcf_alt = revcomp(vcf_alt)

            tmpvcf = tmpdirname / f"{line}.vcf"
            with tmpvcf.open(mode="w") as f_out:
                print(VCF_METALINES, file=f_out)
                header = "\t".join([VCF_HEADER, line])
                print(header, file=f_out)
                print(
                    "\t".join(
                        map(
                            str,
                            [
                                CHROM,
                                pos,
                                MISSING,
                                vcf_ref,
                                vcf_alt,
                                MISSING,
                                MISSING,
                                MISSING,
                                "GT",
                                1,
                            ],
                        )
                    ),
                    file=f_out,
                )

            norm_vcf = tmpdirname.with_suffix(".norm.bcf")
            cmd = f"bcftools norm --check-ref e -f {snakemake.input.reference} -o {norm_vcf} {tmpvcf}"
            args = shlex.split(cmd)
            cp = subprocess.run(args, capture_output=True, text=True)
            if cp.returncode != 0:
                eprint(f"[ERR]: Failed to run bcftools norm for {line}")
                eprint(cp.stderr)

            cmd = f"bcftools index {norm_vcf}"
            args = shlex.split(cmd)
            cp = subprocess.run(args, capture_output=True, text=True)
            if cp.returncode != 0:
                eprint(f"[ERR]: Failed to run bcftools index for {line}")
                eprint(cp.stderr)

            vcfs_to_merge.append(norm_vcf)

    eprint("[INFO]: Merging orphan VCFs...")

    file_list = tmpdirname / "file_list.txt"
    with file_list.open(mode="w") as fp:
        for p in vcfs_to_merge:
            print(str(p), file=fp)

    merged_orphan_vcf = tmpdirname / "orphan.bcf"
    cmd = f"bcftools merge --missing-to-ref --file-list {file_list} -o {merged_orphan_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint("[ERR]: Failed to run bcftools merge for orphan VCF")
        eprint(cp.stderr)

    cmd = f"bcftools index {merged_orphan_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint("[ERR]: Failed to run bcftools index on merged orphan VCF")
        eprint(cp.stderr)

    eprint("[INFO]: Merging orphan VCF with CRyPTIC VCF...")
    cryptic_vcf_file = snakemake.input.cryptic_common_vcf
    merged_vcf = tmpdirname / "merged.bcf"
    cmd = f"bcftools merge --missing-to-ref -o {merged_vcf} {merged_orphan_vcf} {cryptic_vcf_file}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint("[ERR]: Failed to run bcftools merge")
        eprint(cp.stderr)

    cmd = f"bcftools index {merged_orphan_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint("[ERR]: Failed to run bcftools index on merged VCF")
        eprint(cp.stderr)

    eprint("[INFO]: Normalising merged VCF...")
    out_vcf = snakemake.output.vcf
    cmd = f"bcftools norm --check-ref e -f {snakemake.input.reference} -o {out_vcf} {merged_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint(f"[ERR]: Failed to run bcftools norm for merged VCF")
        eprint(cp.stderr)

    cmd = f"bcftools index {out_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint(f"[ERR]: Failed to run bcftools index for output VCF")
        eprint(cp.stderr)

    tmpdir.cleanup()


if __name__ == "__main__":
    main()
