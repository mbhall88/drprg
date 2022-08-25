import sys

sys.stderr = open(snakemake.log[0], "w")

import json
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
    try:
        return data[next(iter(data.keys()))]["susceptibility"]
    except (KeyError, TypeError):
        return data["susceptibility"]


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
        proj = p.parts[-3]
        sample = p.parts[-2]
        run = p.name.split(".")[0]
        tech = snakemake.wildcards.tech

        with open(p) as fp:
            report = load_susceptibility(fp)

        if not report:
            eprint(f"[WARNING] {run} has no susceptibility results")
            continue

        for drug, pred in report.items():
            print(
                DELIM.join(
                    (
                        run,
                        sample,
                        proj,
                        tech,
                        "mykrobe",
                        drug,
                        pred["predict"],
                    )
                ),
                file=fout,
            )
