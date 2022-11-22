import sys

sys.stderr = open(snakemake.log[0], "w")
import matplotlib.pyplot as plt
import seaborn as sns
import upsetplot
import pandas as pd
import numpy as np

plt.style.use("ggplot")


def main():
    samplesheet = snakemake.params.samplesheet

    drop_cols = {"biosample", "bioproject"}
    DRUGS = set(samplesheet.drop(columns=drop_cols).columns)
    pheno_df = (
        samplesheet.drop(columns=drop_cols)
        .melt(ignore_index=False, var_name="drug", value_name="phenotype")
        .reset_index()
    )
    d = {}
    samples_with_pheno = []
    for drug in map(str.lower, DRUGS):
        drug_df = pheno_df.query("drug == @drug").dropna()
        samples = list(drug_df["run"])
        samples_with_pheno.extend(samples)
        if samples:
            d[drug.upper()] = samples

    upset_data = upsetplot.from_contents(d)

    fig, ax = plt.subplots(figsize=(13, 8), dpi=300)
    p = upsetplot.plot(
        upset_data,
        fig=fig,
        sort_by="cardinality",
        orientation="horizontal",
        show_counts=True,
        #     min_subset_size=300
    )
    p["intersections"].set_ylabel("")
    p["intersections"].set_yticks([])
    ax.axis("off")
    p["totals"].set_xticks([])

    for p in snakemake.output.upset_plots:
        fig.savefig(p)

    fig, ax = plt.subplots(dpi=300, figsize=(13, 8), tight_layout=True)
    data = pheno_df.dropna()
    counts = data.groupby(["drug", "phenotype"]).count().reset_index()
    counts["drug"] = [s.capitalize() for s in counts["drug"]]
    sns.barplot(data=counts, x="drug", y="run", hue="phenotype")
    yticks = [
        int(i)
        for i in np.logspace(
            np.log10(max(1, counts["run"].min())), np.log10(counts["run"].max()), num=8
        )
    ]
    ax.set_yscale("log")
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.set(ylabel="Count", xlabel="Drug")
    _ = plt.xticks(rotation=80)

    for p in snakemake.output.barplots:
        fig.savefig(p)


if __name__ == "__main__":
    main()
