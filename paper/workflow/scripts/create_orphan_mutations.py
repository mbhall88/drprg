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
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""

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
    offset = snakemake.params.offset
    reference: dict[str, str] = dict()
    with open(snakemake.input.reference) as fp:
        for line in map(str.rstrip, fp):
            if line and line[0] == ">":
                name = line[1:].split()[0]
                reference[name] = ""
                continue
            reference[name] += line

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
            nucseq = ftr.nucleotide_sequence(reference, offset, offset)
            protseq = ftr.protein_sequence(reference)

            if is_dna:
                mut_start = ((pos + 1 if pos < 1 else pos) - 1) + offset
                mut_end = mut_start + len(ref)
                refseq = nucseq[mut_start:mut_end]
                vcf_ref = refseq
                vcf_pos = mut_start + 1
                vcf_alt = alt
            else:
                mut_start = pos - 1
                refseq = protseq[mut_start]
                dna_start = ((pos * 3 - 2) - 1) + offset
                dna_end = dna_start + 3
                vcf_ref = nucseq[dna_start:dna_end]
                assert CODONTAB[vcf_ref] == ref, line
                vcf_pos = dna_start + 1
                vcf_alt = get_closest_codon(vcf_ref, alt)

            if refseq != ref:
                raise ValueError(
                    f"{line}: ref allele {ref} does not match reference {refseq}"
                )

            tmpvcf = tmpdirname / f"{line}.vcf"
            with tmpvcf.open(mode="w") as f_out:
                print(VCF_METALINES, file=f_out)
                contig_line = f"##contig=<ID={gene},length={len(nucseq)}>"
                print(contig_line, file=f_out)
                header = "\t".join([VCF_HEADER, line])
                print(header, file=f_out)
                print(
                    "\t".join(
                        map(
                            str,
                            [
                                gene,
                                vcf_pos,
                                MISSING,
                                vcf_ref,
                                vcf_alt,
                                MISSING,
                                MISSING,
                                MISSING,
                                "GT",
                                "1/1",
                            ],
                        )
                    ),
                    file=f_out,
                )

            norm_vcf = tmpvcf.with_suffix(".norm.bcf")
            cmd = f"bcftools norm --check-ref e -f {snakemake.input.gene_reference} -o {norm_vcf} {tmpvcf}"
            args = shlex.split(cmd)
            cp = subprocess.run(args, capture_output=True, text=True)
            if cp.returncode != 0:
                eprint(f"[ERR]: Failed to run bcftools norm for {line}")
                eprint(cp.stderr)
                sys.exit(1)

            cmd = f"bcftools index -f {norm_vcf}"
            args = shlex.split(cmd)
            cp = subprocess.run(args, capture_output=True, text=True)
            if cp.returncode != 0:
                eprint(f"[ERR]: Failed to run bcftools index for {line}")
                eprint(cp.stderr)
                sys.exit(1)

            vcfs_to_merge.append(norm_vcf)

    eprint("[INFO]: Merging orphan VCFs...")

    file_list = tmpdirname / "file_list.txt"
    with file_list.open(mode="w") as fp:
        for p in vcfs_to_merge:
            print(str(p), file=fp)

    merged_orphan_vcf = tmpdirname / "orphan.bcf"
    cmd = f"bcftools merge --file-list {file_list} -o {merged_orphan_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint("[ERR]: Failed to run bcftools merge for orphan VCF")
        eprint(cp.stderr)
        sys.exit(1)

    cmd = f"bcftools index -f {merged_orphan_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint("[ERR]: Failed to run bcftools index on merged orphan VCF")
        eprint(cp.stderr)
        sys.exit(1)

    eprint("[INFO]: Normalising merged VCF...")
    norm_vcf = tmpdirname / "norm.bcf"
    cmd = f"bcftools norm --check-ref e -f {snakemake.input.gene_reference} -o {norm_vcf} {merged_orphan_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint(f"[ERR]: Failed to run bcftools norm for merged VCF")
        eprint(cp.stderr)
        sys.exit(1)

    eprint("[INFO]: Sorting normalised VCF...")
    out_vcf = snakemake.output.vcf
    cmd = f"bcftools sort -o {out_vcf} {norm_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint(f"[ERR]: Failed to run bcftools norm for merged VCF")
        eprint(cp.stderr)
        sys.exit(1)

    cmd = f"bcftools index -f {out_vcf}"
    args = shlex.split(cmd)
    cp = subprocess.run(args, capture_output=True, text=True)
    if cp.returncode != 0:
        eprint(f"[ERR]: Failed to run bcftools index for output VCF")
        eprint(cp.stderr)
        sys.exit(1)

    tmpdir.cleanup()


if __name__ == "__main__":
    main()
