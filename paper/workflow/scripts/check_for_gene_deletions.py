"""USAGE: python check_for_gene_deletions.py <indir> <genes.fa> <threads>"""
import sys
from pathlib import Path
from multiprocessing import Pool

N = 100


def missing_genes_in_vcf(path: Path, expected_genes: set[str]):
    run = path.parts[-2]
    genes = set()
    with open(path) as fp:
        for line in map(str.rstrip, fp):
            if not line.startswith("##"):
                break
            if line.startswith("##contig"):
                gene = line.split("=")[-1][:-1]
                genes.add(gene)

    missing_genes = expected_genes - genes
    if missing_genes:
        print(f"{run} is missing: {sorted(missing_genes)}")


def main():
    indir = Path(sys.argv[1])

    fasta_file = sys.argv[2]
    with open(fasta_file) as fp:
        expected_genes = {g[1:].split()[0] for g in fp if g[0] == ">"}
        print(f"Got {len(expected_genes)} expected genes")

    paths = [(p, expected_genes) for p in indir.rglob("pandora_genotyped.vcf")]

    with Pool(int(sys.argv[3])) as pool:
        results = pool.starmap(missing_genes_in_vcf, paths, chunksize=N)


if __name__ == "__main__":
    main()
