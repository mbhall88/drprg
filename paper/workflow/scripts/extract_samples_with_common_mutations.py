import sys
from collections import defaultdict, Counter

import pandas as pd

def eprint(msg: str):
    print(msg, file=sys.stderr)

def load_uids(file: str) -> dict[str, list[str]]:
    d = defaultdict(list)

    seen = set()
    with open(file) as fp:
        for line in map(str.rstrip, fp):
            if not line or line in seen:
                continue
            seen.add(line)
            fields = line.split(".")
            try:
                ix = (fields[5], fields[7])
            except IndexError as err:
                if line.startswith("N"):
                    ix = line
                else:
                    raise err
            d[ix].append(line)

    return d


def main():
    cryptic_samples_file = sys.argv[1]

    cryptic_samples = load_uids(cryptic_samples_file)
    stop = False
    for k, v in cryptic_samples.items():
        if len(v) > 1:
            eprint(f"{k} has {len(v)} uids")
            stop = True

    if stop:
        sys.exit(1)

    popn_samples_file = sys.argv[2]
    popn_samples = load_uids(popn_samples_file)
    for k, v in popn_samples.items():
        if len(v) > 1:
            eprint(f"{k} has {len(v)} uids")
            stop = True

    if stop:
        sys.exit(1)

    mutations_file = sys.argv[3]

    mutations_df = pd.read_csv(mutations_file, low_memory=False)
    mut2id = defaultdict(list)
    counts = Counter()

    for _, row in mutations_df.iterrows():
        mut = row["MUTATION"]
        gene = row["GENE"]
        uid = row["UNIQUEID"]
        mut2id[(gene, mut)].append(uid)
        counts[(gene, mut)] += 1

    min_count = int(sys.argv[4])
    for (gene, mut), count in counts.items():
        if count >= min_count:
            print(gene, mut)

if __name__ == "__main__":
    main()
