import sys

sys.stderr = open(snakemake.log[0], "w")
import json
from pathlib import Path


DELIM = snakemake.params.delim


def load_drprg_report(path: Path) -> dict[str, str]:
    with open(path) as fp:
        data = json.load(fp)
    report = dict()
    for drug, results in data["susceptibility"].items():
        if drug == "NONE":
            continue
        report[drug] = results["predict"]

    return report


def main():
    out_fp = open(snakemake.output.sheet, "w")

    print(
        DELIM.join(["sample", "technology", "w", "k", "drug", "prediction"]),
        file=out_fp,
    )

    for report_path in map(Path, snakemake.input.reports):
        sample = report_path.parts[-2]
        tech = report_path.parts[-3]
        k = report_path.parts[-4][1:]
        w = report_path.parts[-5][1:]

        report = load_drprg_report(report_path)
        for drug, pred in report.items():
            print(DELIM.join([sample, tech, w, k, drug, pred]), file=out_fp)

    out_fp.close()


main()
