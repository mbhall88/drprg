import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use("ggplot")

df = pd.read_csv(snakemake.input.qc)
fig, ax = plt.subplots(figsize=(13, 8), dpi=300, tight_layout=True)
x = "coverage"
sns.histplot(data=df, x=x, ax=ax, log_scale=(True, True))
if snakemake.wildcards.tech == "illumina":
    xticks = [5, 10, 20, 30, 50, 100, 200, 300, 500, 1000]
else:
    xticks = [5, 10, 20, 30, 50, 100, 150, 200]
ax.set_xticklabels(xticks)
ax.set_xticks(xticks)
if snakemake.wildcards.tech == "illumina":
    yticks = [1, 10, 100, 500]
else:
    yticks = [1, 10, 100, 200]
ax.set_yticklabels(yticks)
ax.set_yticks(yticks)
ax.set_xlabel("Read depth")

for p in snakemake.output.plots:
    fig.savefig(p)
