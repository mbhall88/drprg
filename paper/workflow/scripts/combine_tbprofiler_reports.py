import sys

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


def load_susceptibility(stream) -> dict:
    """Extract the susceptibility info from the JSON"""
    data = json.load(stream)
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
            report = load_susceptibility(fp)

        if not report:
            eprint(f"[WARNING] {run} has no susceptibility results")
            continue

        seen_drugs = set()

        for variant in report:
            for info in variant["drugs"]:
                drug = info["drug"]
                if info["confers"] != "resistance" or drug in seen_drugs:
                    continue
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
                        )
                    ),
                    file=fout,
                )
                seen_drugs.add(drug)
