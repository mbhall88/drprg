import sys

sys.stderr = open(snakemake.log[0], "w")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import to_rgba
from statannotations.Annotator import Annotator

plt.style.use("ggplot")
FS = snakemake.params.fontsize
PALETTE = snakemake.params.palette
DPI = snakemake.params.dpi
FIGSIZE = snakemake.params.figsize
STATS_TEST = snakemake.params.stats_test


def legend_without_duplicate_labels(ax, title=""):
    handles, labels = ax.get_legend_handles_labels()
    unique = [
        (h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]
    ]
    ax.legend(*zip(*unique), title=title, fontsize=FS)


def main():
    df = pd.read_csv(snakemake.input.summary)

    # PLOT MEMORY
    mem_fig, mem_ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, tight_layout=True)
    y = "tool"
    hue = "tool"
    x = "max_rss"
    hue_order = sorted(set(df["tool"]))
    pval_fmt = {
        "pvalue_thresholds": [
            [1e-4, "$p\leq0.0001$"],
            [1e-3, "$p\leq0.001$"],
            [1e-2, "$p\leq0.01$"],
            [0.05, "$p\leq0.05$"],
            [1, "ns"],
        ]
    }
    violin_alpha = 0.2
    strip_alpha = 0.5
    orient = "h"

    kwargs = dict(
        data=df,
        x=x,
        y=y,
        hue=hue,
        palette=PALETTE,
        hue_order=hue_order,
        order=hue_order,
        orient=orient,
        dodge=False,
    )

    sns.violinplot(**kwargs, cut=0, inner="quartile", ax=mem_ax)
    for violin in mem_ax.collections:
        violin.set_facecolor(to_rgba(violin.get_facecolor(), alpha=violin_alpha))

    sns.stripplot(**kwargs, alpha=strip_alpha, edgecolor="gray", linewidth=1, ax=mem_ax)

    mem_ax.set_xscale("log")
    ticks = [
        (100, "100MB"),
        (500, "500MB"),
        (1000, "1GB"),
        (2000, "2GB"),
        (3000, "3GB"),
        (4000, "4GB"),
    ]
    mem_ax.set_xticks([t[0] for t in ticks])
    mem_ax.set_xticklabels([t[1] for t in ticks], fontsize=FS)
    mem_ax.set_xlabel("Max. RAM usage", fontsize=FS)
    mem_ax.set_ylabel("")
    mem_ax.tick_params(axis="both", which="major", labelsize=FS)

    pairs = [("drprg", "mykrobe"), ("mykrobe", "tbprofiler"), ("drprg", "tbprofiler")]

    legend_without_duplicate_labels(mem_ax)

    annot = Annotator(mem_ax, pairs, data=df, x=x, y=y, orient=orient, order=hue_order)
    annot.configure(test=STATS_TEST, pvalue_format=pval_fmt)
    annot.apply_test()
    annot.annotate()

    for p in snakemake.output.memory_plots:
        mem_fig.savefig(p)

    # PLOT RUNTIME
    rt_fig, rt_ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, tight_layout=True)
    x = "s"
    kwargs["x"] = x

    sns.violinplot(**kwargs, cut=0, inner="quartile", ax=rt_ax)
    for violin in rt_ax.collections:
        violin.set_facecolor(to_rgba(violin.get_facecolor(), alpha=violin_alpha))

    sns.stripplot(**kwargs, alpha=strip_alpha, edgecolor="gray", linewidth=1, ax=rt_ax)

    rt_ax.set_xscale("log")
    ticks = [
        (60, "1min"),
        (120, "2min"),
        (180, "3min"),
        (300, "5min"),
        (600, "10min"),
        (1800, "30min"),
    ]
    rt_ax.set_xticks([t[0] for t in ticks])
    rt_ax.set_xticklabels([t[1] for t in ticks], fontsize=FS)
    rt_ax.set_xlabel("Runtime", fontsize=FS)
    rt_ax.set_ylabel("")
    rt_ax.tick_params(axis="both", which="major", labelsize=FS)

    legend_without_duplicate_labels(rt_ax)

    annot.reset_configuration()
    annot = Annotator(rt_ax, pairs, data=df, x=x, y=y, orient=orient, order=hue_order)
    annot.configure(test=STATS_TEST, pvalue_format=pval_fmt)
    annot.apply_test()
    annot.annotate()

    for p in snakemake.output.runtime_plots:
        rt_fig.savefig(p)


main()
