import sys
from collections import defaultdict

sys.stderr = open(snakemake.log[0], "w")

import json
import gzip
from pathlib import Path


if str(snakemake.output.report).endswith("csv"):
    DELIM = ","
elif str(snakemake.output.report).endswith("tsv"):
    DELIM = "\t"
else:
    raise NotImplementedError("Don't know output delimiter")


def eprint(msg):
    print(msg, file=sys.stderr)


def load_susceptibility(data: dict) -> dict:
    """Extract the susceptibility info from the JSON"""
    variants = data["dr_variants"]
    return variants


with open(snakemake.output.report, "w") as fout:

    print(
        DELIM.join(
            [
                "run",
                "biosample",
                "bioproject",
                "technology",
                "tool",
                "drug",
                "prediction",
                "lineage",
                "mutations",
            ]
        ),
        file=fout,
    )

    for p in map(Path, snakemake.input.reports):
        proj = p.parts[-5]
        sample = p.parts[-4]
        run = p.parts[-3]
        tech = snakemake.wildcards.tech

        fopen = gzip.open if p.suffix == ".gz" else open

        with fopen(p) as fp:
            data = json.load(fp)
            report = load_susceptibility(data)
            lineage = data["sublin"]

        if not report:
            print(
                DELIM.join(
                    (run, sample, proj, tech, "tbprofiler", "all", "S", lineage, "")
                ),
                file=fout,
            )
            continue

        drug_variants = defaultdict(set)

        for variant in report:
            gene = variant["gene"]
            v = variant["change"]
            mut = f"{gene}_{v}"
            for info in variant["drugs"]:
                drug = info["drug"]
                if info["confers"] != "resistance":
                    continue
                else:
                    drug_variants[drug].add(mut)

        for drug, mutations in drug_variants.items():
            mut_str = ";".join(sorted(mutations))
            print(
                DELIM.join(
                    (
                        run,
                        sample,
                        proj,
                        tech,
                        "tbprofiler",
                        drug,
                        "R",
                        lineage,
                        mut_str,
                    )
                ),
                file=fout,
            )
