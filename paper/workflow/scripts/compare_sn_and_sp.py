import sys

sys.stderr = open(snakemake.log[0], "w")

from collections import defaultdict, Counter
from dataclasses import dataclass
from enum import Enum
from itertools import product
from math import sqrt
from typing import Tuple, Optional

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from scipy import stats

plt.style.use("ggplot")
FIGSIZE = snakemake.params.figsize
DPI = snakemake.params.dpi
CMAP = plt.get_cmap("Set2").colors


class Prediction(Enum):
    Resistant = "R"
    Susceptible = "S"
    MinorResistance = "r"
    Unknown = "U"
    Failed = "F"

    def __str__(self) -> str:
        return self.value


class Classification(Enum):
    TruePositive = "TP"
    FalsePositive = "FP"
    TrueNegative = "TN"
    FalseNegative = "FN"

    def __str__(self) -> str:
        return self.value


class Classifier:
    def __init__(
        self,
        minor_is_susceptible: bool = False,
        unknown_is_resistant: bool = False,
        failed_is_resistant: bool = False,
    ):
        self.minor_is_susceptible = minor_is_susceptible
        self.unknown_is_resistant = unknown_is_resistant
        self.failed_is_resistant = failed_is_resistant
        self.susceptible = {Prediction.Susceptible}
        self.resistant = {Prediction.Resistant}
        if self.minor_is_susceptible:
            self.susceptible.add(Prediction.MinorResistance)
        else:
            self.resistant.add(Prediction.MinorResistance)

        if self.unknown_is_resistant:
            self.resistant.add(Prediction.Unknown)
        else:
            self.susceptible.add(Prediction.Unknown)

        if self.failed_is_resistant:
            self.resistant.add(Prediction.Failed)
        else:
            self.susceptible.add(Prediction.Failed)

    def from_predictions(
        self, y_true: Prediction, y_pred: Prediction
    ) -> Classification:
        if y_true in self.susceptible:
            expected_susceptible = True
        elif y_true in self.resistant:
            expected_susceptible = False
        else:
            raise NotImplementedError(f"Don't know how to classify {y_true} calls yet")

        if y_pred in self.susceptible:
            called_susceptible = True
        elif y_pred in self.resistant:
            called_susceptible = False
        else:
            raise NotImplementedError(f"Don't know how to classify {y_pred} calls yet")

        if expected_susceptible and called_susceptible:
            return Classification.TrueNegative
        elif expected_susceptible and not called_susceptible:
            return Classification.FalsePositive
        elif not expected_susceptible and not called_susceptible:
            return Classification.TruePositive
        else:
            return Classification.FalseNegative


@dataclass
class ConfusionMatrix:
    tp: int = 0
    tn: int = 0
    fp: int = 0
    fn: int = 0

    def ravel(self) -> Tuple[int, int, int, int]:
        """Return the matrix as a flattened tuple.
        The order of return is TN, FP, FN, TP
        """
        return self.tn, self.fp, self.fn, self.tp

    def as_matrix(self) -> np.ndarray:
        """Returns a 2x2 matrix [[TN, FP], [FN, TP]]"""
        return np.array([[self.tn, self.fp], [self.fn, self.tp]])

    def num_positive(self) -> int:
        """Number of TPs and FNs - i.e. actual condition positive"""
        return self.tp + self.fn

    def num_negative(self) -> int:
        """Number of TNs and FPs - i.e. actual condition negative"""
        return self.tn + self.fp

    def sensitivity(self) -> Tuple[float, float, float]:
        """Also known as recall and true positive rate (TPR)"""
        try:
            sn = self.tp / (self.tp + self.fn)
            lwr_bound, upr_bound = confidence_interval(n_s=self.tp, n_f=self.fn)
            return sn, lwr_bound, upr_bound
        except ZeroDivisionError:
            return None, None, None

    def specificity(self) -> Tuple[float, float, float]:
        """Also known as selectivity and true negative rate (TNR)"""
        try:
            sp = self.tn / (self.tn + self.fp)
            lwr_bound, upr_bound = confidence_interval(n_s=self.tn, n_f=self.fp)
            return sp, lwr_bound, upr_bound
        except ZeroDivisionError:
            return None, None, None

    def mcc(self) -> Optional[float]:
        """Matthews correlation coefficient"""
        numerator = (self.tp * self.tn) - (self.fp * self.fn)
        denominator = (
            (self.tp + self.fp)
            * (self.tp + self.fn)
            * (self.tn + self.fp)
            * (self.tn + self.fn)
        )
        try:
            return numerator / sqrt(denominator)
        except ZeroDivisionError:
            return None

    @staticmethod
    def from_series(s: pd.Series) -> "ConfusionMatrix":
        tp = s.get("TP", 0)
        fp = s.get("FP", 0)
        fn = s.get("FN", 0)
        tn = s.get("TN", 0)
        return ConfusionMatrix(tp=tp, fn=fn, fp=fp, tn=tn)


def confidence_interval(n_s: int, n_f: int, conf: float = 0.95) -> Tuple[float, float]:
    """Calculate the Wilson score interval.
    Equation take from https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval
    n_s: Number of successes or, in the case of confusion matrix statistics, the numerator
    n_f: Number of failures or, in the case of confusion matrix statistics, the denominator minus the numerator
    conf: the confidence level. i.e. 0.95 is 95% confidence
    """
    n = n_f + n_s
    z = stats.norm.ppf(1 - (1 - conf) / 2)  # two-sided
    z2 = z**2
    nz2 = n + z2
    A = (n_s + (0.5 * z2)) / nz2
    B = z / nz2
    C = sqrt(((n_s * n_f) / n) + (z2 / 4))
    CI = B * C
    return A - CI, A + CI


def round_up_to_base(x, base=10):
    return int(x + (base - x) % base)


def round_down_to_base(x, base=10):
    return int(x - (x % base))


def main():
    ignore_drugs = snakemake.params.ignore_drugs
    frames = []
    for p in snakemake.input.summary_files:
        frames.append(pd.read_csv(p))
    calls = pd.concat(frames)
    calls["drug"] = calls.drug.str.lower()
    tools = sorted(set(calls["tool"]))
    runs = set(calls["run"])
    calls.query("drug not in @ignore_drugs", inplace=True)
    calls.set_index(
        ["run", "tool", "drug"], verify_integrity=True, inplace=True, drop=False
    )

    phenotypes = pd.read_csv(
        snakemake.input.phenotypes, index_col="run", low_memory=False
    )
    minor_is_susceptible = snakemake.params.minor_is_susceptible
    unknown_is_resistant = snakemake.params.unknown_is_resistant
    failed_is_resistant = snakemake.params.failed_is_resistant
    classifier = Classifier(
        unknown_is_resistant=unknown_is_resistant,
        minor_is_susceptible=minor_is_susceptible,
        failed_is_resistant=failed_is_resistant,
    )

    drugs = set()
    pheno_clf = []
    for drug, tool, run in product(set(calls["drug"]), tools, runs):
        try:
            ph = phenotypes.at[run, drug]
            if pd.isna(ph):
                continue
            else:
                truth = Prediction(ph)
        except KeyError:
            print(f"[WARN]: {run} has no phenotype column for {drug}")
            continue

        ix = (run, tool, drug)
        drugs.add(drug)
        try:
            pred = calls.at[ix, "prediction"]
        except KeyError:
            # tbprofiler doesnt explicitly report S
            assert tool == "tbprofiler", ix
            pred = "S"
        clf = classifier.from_predictions(truth, Prediction(pred))

        pheno_clf.append((run, drug, str(clf), tool))
        cols = ["run", "drug", "classification", "tool"]
    df = pd.DataFrame(pheno_clf, columns=cols)

    df.to_csv(snakemake.output.classification, index=False)

    cms = defaultdict()

    low_pheno_drugs = set()

    for drug, tool in product(drugs, tools):
        s = df.query("drug == @drug and tool == @tool").value_counts(
            subset=["classification"]
        )
        cm = ConfusionMatrix.from_series(s)
        cms[(drug, tool)] = cm
        n = sum(cm.ravel())
        if n < snakemake.params.min_num_phenotypes:
            low_pheno_drugs.add(drug)

    metrics = []
    for (drug, tool), cm in cms.items():
        sn = cm.sensitivity()[0]
        sp = cm.specificity()[0]
        metrics.append((drug, tool, sn, sp))

    summary_cols = [
        "drug",
        "tool",
        "Sensitivity",
        "Specificity",
    ]

    summary = pd.DataFrame(
        metrics,
        columns=summary_cols,
    ).melt(id_vars=["drug", "tool"], var_name="metric")

    table = (
        summary.set_index(["drug", "tool", "metric"])["value"].unstack().reset_index()
    )

    for clf in ["TP", "FP", "FN", "TN"]:
        table[clf] = 0

    for i, row in table.iterrows():
        ix = (row["drug"], row["tool"])
        cm = cms[ix]
        table.at[i, "TP"] = cm.tp
        table.at[i, "FP"] = cm.fp
        table.at[i, "TN"] = cm.tn
        table.at[i, "FN"] = cm.fn

    for k in ["drug", "tool"]:
        table[k] = table[k].str.capitalize()

    table.fillna("-", inplace=True)
    summary_cols = [
        "drug",
        "tool",
        "Sensitivity",
        "Specificity",
        "TP",
        "TN",
        "FN",
        "FP",
    ]
    table = table[summary_cols]

    rows = []
    ci_str = (
        lambda tup: f"{tup[0]:.1%} ({tup[1] * 100:.1f}-{tup[2]:.1%})"
        if tup[0] is not None
        else "-"
    )
    for i, row in table.iterrows():
        cm = ConfusionMatrix(tp=row["TP"], fp=row["FP"], tn=row["TN"], fn=row["FN"])
        sn = cm.sensitivity()
        sp = cm.specificity()
        mcc = cm.mcc()
        if mcc is not None:
            mcc = round(mcc, 3)
        else:
            mcc = "-"
        fn_str = f"{cm.fn}({cm.num_positive()})"
        fp_str = f"{cm.fp}({cm.num_negative()})"
        rows.append(
            (
                row["drug"].capitalize(),
                row["tool"],
                fn_str,
                fp_str,
                ci_str(sn),
                ci_str(sp),
                mcc,
            )
        )
    pretty_cols = [
        "Drug",
        "Tool",
        "FN(R)",
        "FP(S)",
        "Sensitivity (95% CI)",
        "Specificity (95% CI)",
        "MCC",
    ]
    table = pd.DataFrame(rows, columns=pretty_cols)

    table.to_csv(snakemake.output.table, index=False)

    sn_data = []
    sp_data = []
    for drug, tool in product(drugs, tools):
        if drug in low_pheno_drugs:
            continue

        s = df.query("drug == @drug and tool == @tool").value_counts(
            subset=["classification"]
        )
        cm = ConfusionMatrix.from_series(s)
        sn = cm.sensitivity()
        sp = cm.specificity()
        sn_data.append((drug, tool, *sn))
        sp_data.append((drug, tool, *sp))

    sn_df = pd.DataFrame(sn_data, columns=["drug", "tool", "value", "lower", "upper"])
    sp_df = pd.DataFrame(sp_data, columns=["drug", "tool", "value", "lower", "upper"])

    s = """AMK amikacin
    CAP capreomycin
    CFX ciprofloxacin
    DLM delamanid
    EMB ethambutol
    ETO ethionamide
    INH isoniazid
    KAN kanamycin
    LFX levofloxacin
    LZD linezolid
    MFX moxifloxacin
    OFX ofloxacin
    PZA pyrazinamide
    RIF rifampicin
    STM streptomycin"""
    long2short = dict()
    for line in s.splitlines():
        ab, d = line.split()
        long2short[d] = ab

    short2long = dict()
    for line in s.splitlines():
        ab, d = line.split()
        short2long[ab] = d

    first_line = [short2long[d] for d in ["INH", "RIF", "EMB", "PZA"]]
    fluoroquinolones = [short2long[d] for d in ["LFX", "MFX", "OFX", "CFX"]]  # group A
    macrolides = [
        short2long[d] for d in ["AMK", "CAP", "KAN", "STM"]
    ]  # group B - second-line injectables
    other = [short2long[d] for d in ["ETO", "LZD", "DLM"]]  # group C and D
    drug_order = [
        d
        for d in [*first_line, *fluoroquinolones, *macrolides, *other]
        if d in drugs and d not in low_pheno_drugs
    ]

    def sort_drugs(a):
        xs = drug_order
        out = []
        c = Counter()
        for x in a:
            i = xs.index(x)
            d = xs[i]
            c[d] += 1
            out.append((i, c[d]))
        return out

    sn_df = sn_df.sort_values(by="drug", key=sort_drugs).reset_index(drop=True)
    sp_df = sp_df.sort_values(by="drug", key=sort_drugs).reset_index(drop=True)

    # PLOT
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)

    # plot details
    bar_width = 0.2
    epsilon = 0.05
    fontsize = 12
    rotate = 0
    dodge = 0.1
    capsize = 2
    marker_alpha = 0.8
    marker_size = 7
    edge_alpha = 0.7
    edgecol = to_rgba("black", alpha=edge_alpha)
    edge_line_width = 1
    sn_marker = snakemake.params.sn_marker
    sp_marker = snakemake.params.sp_marker
    legend_marker_size = 8

    i = -1
    all_positions = []
    leghandles = []

    for tool in tools:
        i += 1
        positions = [
            (p - 1) + ((bar_width + epsilon) * i) for p in np.arange(len(drugs))
        ]

        all_positions.append(positions)

        colour = to_rgba(CMAP[i], alpha=marker_alpha)
        label = tool

        tool_sp_df = sp_df.query("tool==@tool")
        sp_ys = tool_sp_df["value"] * 100
        sp_lb = sp_ys - tool_sp_df["lower"] * 100
        sp_ub = tool_sp_df["upper"] * 100 - sp_ys
        sp_ub = [min(100, x) for x in sp_ub]

        plotprops = dict(
            label=label,
            color=colour,
            capsize=capsize,
            elinewidth=edge_line_width,
            mec=edgecol,
            markersize=marker_size,
        )

        sp_bar = ax.errorbar(
            x=positions, y=sp_ys, yerr=[sp_lb, sp_ub], fmt=sp_marker, **plotprops
        )

        tool_sn_df = sn_df.query("tool==@tool")
        sn_ys = tool_sn_df["value"] * 100
        sn_lb = sn_ys - tool_sn_df["lower"] * 100
        sn_ub = tool_sn_df["upper"] * 100 - sn_ys

        sn_bar = ax.errorbar(
            x=[p - dodge for p in positions],
            y=sn_ys,
            yerr=[sn_lb, sn_ub],
            fmt=sn_marker,
            **plotprops,
        )

        h = mlines.Line2D(
            [],
            [],
            color=CMAP[i],
            marker=sp_marker,
            markersize=legend_marker_size,
            label=label,
        )
        leghandles.append(h)

    labels = [long2short[d] for d in drug_order]
    label_pos = [np.mean(ps) for ps in zip(*all_positions)]
    plt.xticks(label_pos, labels, rotation=rotate, fontsize=fontsize)
    ax.set_ylabel("Sensitivity/Specificity (%)")
    yticks = [0, 10, 30, 40, 50, 60, 70, 80, 85, 90, 95, 100]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.set_xticks(label_pos)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=rotate, fontsize=fontsize)
    ax.tick_params("both", labelsize=fontsize)

    leghandles.append(
        mlines.Line2D(
            [],
            [],
            color="black",
            marker=sn_marker,
            markersize=legend_marker_size,
            label="Sensitivity",
        )
    )
    leghandles.append(
        mlines.Line2D(
            [],
            [],
            color="black",
            marker=sp_marker,
            markersize=legend_marker_size,
            label="Specificity",
        )
    )

    legend_props = dict(
        loc="upper center",
        prop=dict(size=fontsize),
        frameon=False,
        ncol=len(leghandles),
    )
    ax.legend(handles=leghandles, bbox_to_anchor=(0.5, 1.05), **legend_props)

    ax.grid(
        False,
        axis="x",
    )
    # draw line between drug bars
    for xpos in all_positions[-1]:
        vpos = xpos + epsilon + bar_width
        ax.axvline(vpos, color="white", linestyle="-", alpha=1, linewidth=6)

    xlim = (all_positions[0][0] - (epsilon + bar_width), vpos + epsilon)
    ax.set_xlim(xlim)
    _ = ax.set_xlabel("Drug", fontsize=fontsize + 2)
    plt.tight_layout()

    for plot in snakemake.output.plots:
        fig.savefig(plot)


if __name__ == "__main__":
    main()
