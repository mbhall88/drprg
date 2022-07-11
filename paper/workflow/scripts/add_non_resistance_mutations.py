import re
import sys
from pathlib import Path
from typing import Dict, Tuple

sys.stderr = open(snakemake.log[0], "w")

TRANSLATE = str.maketrans("ATGC", "TACG")


def eprint(msg: str):
    print(msg, file=sys.stderr)


def load_strands(path: Path) -> Dict[str, str]:
    strands = dict()
    with open(path) as fp:
        for line in map(str.rstrip, fp):
            if line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            ftr = fields[2]
            if ftr != "gene":
                continue
            extras = fields[8].split(";")
            gene = ""
            for e in extras:
                if e.startswith("Name=") or e.startswith("gene="):
                    gene = e.split("=")[-1]
                    break
            if not gene:
                raise ValueError(f"Couldn't find the gene name for row {line}")
            strands[gene] = fields[6]

    return strands


def revcomp(s: str) -> str:
    return complement(s)[::-1]


def complement(s: str) -> str:
    return s.upper().translate(TRANSLATE)


def split_var_name(name: str) -> Tuple[str, int, str]:
    items = re.match(r"([A-Z]+)([-0-9]+)([A-Z/\*]+)", name, re.I).groups()
    return items[0], int(items[1]), items[2]


def main():
    gene_strands = load_strands(snakemake.input.features)
    exclude = set(snakemake.params.exclude)
    out_fp = open(snakemake.output.panel, "w")

    panel_genes: set[str] = set()
    with open(snakemake.input.panel) as fp:
        for line in fp:
            gene = line.split("\t")[0]
            panel_genes.add(gene)
            out_fp.write(line)

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
                if gene_strands[gene] == "-" and alpha == "DNA" and not is_promotor_mut:
                    ref, pos, alt = split_var_name(var)
                    ref = revcomp(ref)
                    alt = revcomp(alt)
                    pos -= len(ref) - 1
                    var = f"{ref}{pos}{alt}"

                print(
                    "\t".join([gene, var, alpha, snakemake.params.no_drug]),
                    file=out_fp,
                )
                c += 1

    eprint(f"Added {c} non-resistant mutations to panel")

    out_fp.close()


main()
