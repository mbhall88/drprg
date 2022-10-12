import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path

DELIM = snakemake.params.delim


def main():
    header_written = False

    with open(snakemake.output.summary, "w") as fp_out:
        header = ["techology", "tool", "run"]

        for file in map(Path, snakemake.input.bench):
            run = file.stem
            tool = file.parts[-5]

            with open(file) as fp_in:
                first_line = next(fp_in).rstrip()
                if not header_written:
                    columns = first_line.split("\t")
                    header.extend(columns)
                    print(DELIM.join(header), file=fp_out)
                    header_written = True

                for line in map(str.rstrip, fp_in):
                    fields = line.split("\t")
                    row = [snakemake.wildcards.tech, tool, run, *fields]
                    print(DELIM.join(row), file=fp_out)


main()
