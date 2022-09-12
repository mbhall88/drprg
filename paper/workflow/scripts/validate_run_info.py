import sys

sys.stderr = open(snakemake.log[0], "w")

import json
from pathlib import Path

DELIM = snakemake.params.delim
NL = "\n"
PAIRED = "PAIRED"
SINGLE = "SINGLE"


class PlatformError(Exception):
    pass


class LayoutError(Exception):
    pass


def eprint(msg):
    print(msg, file=sys.stderr)


def main():
    out_data = []

    p = Path(snakemake.input.info)
    d = p.parent
    run = d.parts[-1]
    with open(p) as fp:
        data = json.load(fp)[0]

        platform = data["instrument_platform"]
        tech = snakemake.wildcards.tech.upper()
        if tech not in platform:
            raise PlatformError(
                f"Expected {tech}, but got platform {platform} for run {run}"
            )

        layout = data["library_layout"]
        model = data["instrument_model"]
        fastqs = list(d.glob("*.fastq.gz"))

        if not all(p.exists() for p in fastqs):
            raise FileNotFoundError(f"One or more fastqs don't exist {fastqs}")

        if layout == PAIRED and len(fastqs) < 2:
            raise LayoutError(f"Expected paired fastqs, but only got {fastqs}")

        if tech == "nanopore":
            assert layout == SINGLE
            if len(fastqs) > 1:
                raise NotImplementedError(f"Got more than one fastq for Nanopore: {d}")

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

    out_data.append((run, files, model))

    with open(snakemake.output.run_info, "w") as fp:
        print(f"run{DELIM}layout{DELIM}model", file=fp)
        for t in out_data:
            print(DELIM.join(t), file=fp)


if __name__ == "__main__":
    main()
