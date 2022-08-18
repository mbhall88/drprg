import sys

sys.stderr = open(snakemake.log[0], "w")

import json
from pathlib import Path

DELIM = snakemake.params.delim
NL = "\n"
LAYOUT = {"SINGLE": 1, "PAIRED": 2}


def eprint(msg):
    print(msg, file=sys.stderr)


def main():
    out_data = []

    wrong_platform = []

    for d in map(Path, snakemake.input.dirs):
        run = d.parts[-1]
        p = d / "fastq-run-info.json"
        with open(p) as fp:
            data = json.load(fp)[0]

            platform = data["instrument_platform"]
            tech = snakemake.wildcards.tech.upper()
            if tech not in platform:
                wrong_platform.append(run)
                eprint(f"Expected {tech}, but got platform {platform} for run {run}")
                continue

            layout = LAYOUT.get(data["library_layout"])
            fastqs = [d / Path(p).name for p in data["fastq_ftp"].split(";")]
            if layout is None or layout != len(fastqs):
                raise ValueError(
                    f"Library layout ({data['library_layout']}) does not match the number of fastq files ({len(fastqs)}) for accession {d}"
                )

            if not all(p.exists() for p in fastqs):
                raise FileNotFoundError(f"One or more fastqs don't exist {fastqs}")

            files = ";".join(map(str, fastqs))

            if data["tax_id"] != "1773":
                eprint(f"[WARNING]: Got non-MTB tax ID for {run} - {data['tax_id']}")

        out_data.append((run, files))

    if wrong_platform:
        raise ValueError(
            f"Got an unexpected platform for the following run accessions:\n{NL.join(wrong_platform)}"
        )

    with open(snakemake.output.run_info, "w") as fp:
        print(f"run{DELIM}layout", file=fp)
        for t in out_data:
            print(DELIM.join(t), file=fp)


if __name__ == "__main__":
    main()
