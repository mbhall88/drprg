import sys
from itertools import chain

sys.stderr = open(snakemake.log[0], "w")

from collections import defaultdict

import pandas as pd


def eprint(msg: str):
    print(msg, file=sys.stderr)


def load_uids(file: str) -> dict[tuple[str, str], list[str]]:
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


def flatten(xs):
    """Flatten one level of nesting"""
    return chain.from_iterable(xs)


def main():
    with open(snakemake.input.mutations) as fp:
        mutations = {m for m in map(str.rstrip, fp) if m}

    cryptic_samples_file = snakemake.input.cryptic_samples

    cryptic_samples = load_uids(cryptic_samples_file)

    popn_samples_file = snakemake.input.popn_samples
    popn_samples = load_uids(popn_samples_file)

    mutations_file = snakemake.input.mutations2samples

    mutations_df = pd.read_csv(mutations_file, low_memory=False)
    muts = set(mutations_df["MUTATION"])
    c = 0
    mut2ids = defaultdict(set)

    for v in mutations:
        gene = v.split("_")[0]
        assert "(" not in v, v

        vnt = v.split("_", maxsplit=1)[1]
        if vnt in muts:
            c += 1
            frame = mutations_df.query("MUTATION==@vnt and GENE==@gene")
            mut2ids[v].update(set(frame["UNIQUEID"]))
            continue
        # x_vnt = vnt[:-1] + "X"
        # if x_vnt in muts:
        #     c += 1
        #     frame = mutations_df.query("MUTATION==@x_vnt and GENE==@gene")
        #     mut2ids[v].update(set(frame["UNIQUEID"]))
        #     continue
        # o_vnt = vnt[:-1] + "O"
        # if o_vnt in muts:
        #     c += 1
        #     frame = mutations_df.query("MUTATION==@o_vnt and GENE==@gene")
        #     mut2ids[v].update(set(frame["UNIQUEID"]))
        #     continue
        if "del" in vnt or "ins" in vnt:
            pos = vnt.split("_")[0]
            indel_v = f"{pos}_indel"
            if indel_v in muts:
                c += 1
                frame = mutations_df.query("MUTATION==@indel_v and GENE==@gene")
                mut2ids[v].update(set(frame["UNIQUEID"]))
                continue
        mut2ids[v].update([])

    eprint(f"{c}/{len(mutations)} requested mutations are in the CRyPTIC mutations")

    sample2muts = defaultdict(list)
    no_sample = set()
    for v, ids in mut2ids.items():
        for uid in ids:
            key = (uid.split(".")[3], uid.split(".")[5])
            cryptic = cryptic_samples.get(key, cryptic_samples.get((key[1], key[0])))
            if cryptic is None:
                no_sample.add(uid)
                continue
            sample2muts[cryptic[0]].append(v)

    muts_covered = set(flatten(sample2muts.values()))
    n_muts_covered = len(muts_covered)
    eprint(
        f"{n_muts_covered}/{len(mutations)} requested mutations covered by {len(sample2muts)} samples"
    )

    samples_to_use = set()
    muts_seen = set()
    for sample in sorted(sample2muts, key=lambda k: len(sample2muts[k]), reverse=True):
        for v in sample2muts[sample]:
            if v in muts_seen:
                continue
            muts_seen.add(v)
            samples_to_use.add(sample)
            if len(muts_seen) == len(muts_covered):
                break

    samples_already_used_in_popn_vcf = set(flatten(popn_samples.values()))
    samples_to_use -= samples_already_used_in_popn_vcf
    eprint(f"{len(samples_to_use)} additional samples to be used with requested variants")

    with open(snakemake.output.samples, "w") as fp:
        for s in samples_to_use:
            print(s, file=fp)

    mutations_with_no_cyptic_sample = mutations - muts_covered

    eprint(f"{len(mutations_with_no_cyptic_sample)} mutations with no cryptic sample")

    with open(snakemake.output.orphan_mutations, "w") as fp:
        for m in mutations_with_no_cyptic_sample:
            print(m, file=fp)


if __name__ == "__main__":
    main()
