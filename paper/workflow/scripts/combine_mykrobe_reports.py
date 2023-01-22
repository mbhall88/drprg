import sys
from itertools import chain

sys.stderr = open(snakemake.log[0], "w")

import json
import gzip
from pathlib import Path
from multiprocessing import Pool

N = 100  # number of files to pass to a process

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


def load_report(p: Path) -> list[str]:
    proj = p.parts[-3]
    sample = p.parts[-2]
    run = p.name.split(".")[0]
    tech = snakemake.wildcards.tech

    fopen = gzip.open if p.suffix == ".gz" else open

    with fopen(p) as fp:
        report = load_susceptibility(fp)

    if not report:
        raise ValueError(f"{run} has no susceptibility results")

    rows = []

    for drug, pred in report.items():
        evidence = pred.get("called_by", dict())
        row = DELIM.join(
            (
                run,
                sample,
                proj,
                tech,
                "mykrobe",
                drug,
                pred["predict"],
                evidence_to_str(evidence),
            )
        )
        rows.append(row)

    return rows


def evidence_to_str(evidence: dict[str, dict]) -> str:
    mut_strs = []
    for variant in evidence:
        mut_strs.append(variant.split("-")[0])

    return ";".join(mut_strs)


def flatten(list_of_lists):
    """Flatten one level of nesting"""
    return chain.from_iterable(list_of_lists)


def main():
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
                    "mutations",
                ]
            ),
            file=fout,
        )

        mykrobe_reports = list(map(Path, snakemake.input.reports))

        with Pool(snakemake.threads) as pool:
            rows = pool.map(load_report, mykrobe_reports, chunksize=N)

        for row in flatten(rows):
            print(row, file=fout)


if __name__ == "__main__":
    main()
