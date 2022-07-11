import sys

sys.stderr = open(snakemake.log[0], "w")


def eprint(msg: str):
    print(msg, file=sys.stderr)


def main():
    exclude = set(snakemake.params.exclude)
    out_fp = open(snakemake.output.panel, "w")

    panel_genes: set[str] = set()
    with open(snakemake.input.panel) as fp:
        for line in fp:
            gene = line.split("\t")[0]
            panel_genes.add(gene)
            out_fp.write(line)

    eprint(f"Loaded {len(panel_genes)} genes from panel")

    c = 0
    with open(snakemake.input.known) as fp:
        for line in fp:
            fields = line.rstrip().split("\t")
            gene = fields[0]
            var = fields[1]
            mut = f"{gene}_{var}"
            if gene not in panel_genes or mut in exclude:
                continue

            grade = int(fields[-1])
            if grade in snakemake.params.keep_grades:
                print(
                    "\t".join([gene, fields[1], fields[2], snakemake.params.no_drug]),
                    file=out_fp,
                )
                c += 1

    eprint(f"Added {c} non-resistant mutations to panel")

    out_fp.close()


main()
