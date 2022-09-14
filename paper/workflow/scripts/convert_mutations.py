import sys

sys.stderr = open(snakemake.log[0], "w")
import json

DELIM = "\t"


def main():
    with open(snakemake.input.var2drug) as fp:
        var2drug = json.load(fp)

    with open(snakemake.input.panel) as in_fp, open(
        snakemake.output.panel, "w"
    ) as out_fp:
        for line in map(str.rstrip, in_fp):
            gene, var, alphabet = line.split(DELIM)
            mut = f"{gene}_{var}"
            drugs: list[str] = var2drug[mut]

            if mut == snakemake.params["old"]:
                mut = snakemake.params["new"]
                gene, var = mut.split("_")
                alphabet = snakemake.params["new_alphabet"]

            print(DELIM.join([gene, var, alphabet, ",".join(drugs)]), file=out_fp)


main()
