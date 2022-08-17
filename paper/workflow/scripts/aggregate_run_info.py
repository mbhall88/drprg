import sys

sys.stderr = open(snakemake.log[0], "w")

import json
from pathlib import Path

DELIM = snakemake.params.delim


def eprint(msg):
    print(msg, file=sys.stderr)


out_data = []

for d in map(Path, snakemake.input.dirs):
    run = d.parts[-1]
    p = d / "fastq-run-info.json"
    with open(p) as fp:
        data = json.load(fp)[0]

        platform = data["instrument_platform"]
        tech = snakemake.wildcards.tech.upper()
        if tech not in platform:
            raise ValueError(f"Expected {tech}, but got platform {platform}")

        layout = data["library_layout"]
        r1 = d / f"{run}_1.fastq.gz"
        fq = d / f"{run}.fastq.gz"
        if r1.exists():
            r2 = d / f"{run}_2.fastq.gz"
            if not r2.exists():
                raise FileNotFoundError(f"Found R1 but not R2 file {r2}")
            if layout == "SINGLE":
                eprint(
                    f"[WARNING]: Found PE reads for {run} but run info says SE library...using PE..."
                )
            if fq.exists():
                eprint(f"[INFO]: Removing orphan SE file for {run} PE accession {fq}")
            files = f"{r1};{r2}"
        elif not fq.exists():
            raise FileNotFoundError(f"Expected single-end library {fq}")
        else:  # no PE reads, but SE reads found
            if layout == "PAIRED":
                eprint(
                    f"[WARNING]: Found SE reads for {run} but run info says PE library...using SE..."
                )
            files = str(fq)

        if data["tax_id"] != "1773":
            eprint(f"[WARNING]: Got non-MTB tax ID for {run} - {data['tax_id']}")

    out_data.append((run, files))

with open(snakemake.output.run_info, "w") as fp:
    print(f"run{DELIM}layout", file=fp)
    for t in out_data:
        print(DELIM.join(t), file=fp)
