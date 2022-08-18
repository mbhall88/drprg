import sys

sys.stderr = open(snakemake.log[0], "w")

import json
import os
from pathlib import Path

DELIM = snakemake.params.delim
NL = "\n"
PAIRED = "PAIRED"
SINGLE = "SINGLE"


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

            layout = data["library_layout"]
            fastqs = [d / Path(p).name for p in data["fastq_ftp"].split(";")]

            if not all(p.exists() for p in fastqs):
                raise FileNotFoundError(f"One or more fastqs don't exist {fastqs}")

            if layout == PAIRED and len(fastqs) < 2:
                raise FileNotFoundError(
                    f"Expected paired fastqs, but only got {fastqs}"
                )

            if tech == "nanopore":
                assert layout == SINGLE
                if len(fastqs) > 1:
                    raise NotImplementedError(
                        f"Got more than one fastq for Nanopore: {d}"
                    )
                fq = d / f"{run}.fastq.gz"
                if fq != fastqs[0]:
                    fastqs[0].rename(fq)

                files = str(fq)
            else:
                loner = d / f"{run}.fastq.gz"
                r1 = d / f"{run}_1.fastq.gz"
                r2 = d / f"{run}_2.fastq.gz"
                if layout == PAIRED:
                    assert r1.exists(), d
                    assert r2.exists(), d
                    files = f"{r1};{r2}"
                    if loner.exists():
                        loner.unlink()
                else:
                    if r1.exists():
                        files = str(r1)
                        if loner.exists():
                            loner.unlink()
                    elif loner.exists():
                        files = str(loner)
                    else:
                        raise FileNotFoundError(
                            f"Couldn't find an expected single-end file for {d}"
                        )

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
