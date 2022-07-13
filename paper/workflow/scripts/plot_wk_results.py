import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from enum import Enum
from dataclasses import dataclass
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from matplotlib import lines
from matplotlib.lines import Line2D
import seaborn as sns


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

    def ravel(self) -> tuple[int, int, int, int]:
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

    def sensitivity(self) -> float | None:
        """Also known as recall and true positive rate (TPR)"""
        try:
            return self.tp / self.num_positive()
        except ZeroDivisionError:
            return None

    def specificity(self) -> float | None:
        """Also known as selectivity and true negative rate (TNR)"""
        try:
            return self.tn / self.num_negative()
        except ZeroDivisionError:
            return None

    @staticmethod
    def from_series(s: pd.Series) -> "ConfusionMatrix":
        tp = s.get("TP", 0)
        fp = s.get("FP", 0)
        fn = s.get("FN", 0)
        tn = s.get("TN", 0)
        return ConfusionMatrix(tp=tp, fn=fn, fp=fp, tn=tn)


def main():
    df = pd.read_csv(snakemake.input.sheet)
    df["drug"] = df["drug"].apply(str.capitalize)
    valid_samples = set(df["sample"])

    phenotypes = (
        pd.read_csv(snakemake.input.phenotypes)
        .melt(id_vars=["sample"], var_name="drug", value_name="phenotype")
        .query("sample in @valid_samples")
    )
    phenotypes["drug"] = phenotypes["drug"].apply(str.capitalize)
    # drop lpa phenotypes
    phenotypes = phenotypes[~phenotypes["drug"].str.contains("-lpa")]
    # remove non R/S phenotypes
    arr = []
    for r in phenotypes["phenotype"]:
        if pd.isna(r):
            arr.append(r)
        elif r.upper() in ("R", "S"):
            arr.append(r.upper())
        else:
            arr.append(None)
    phenotypes["phenotype"] = arr
    phenotypes.set_index(["sample", "drug"], drop=False, inplace=True)
    phenotypes.sort_index(inplace=True)
    # drop rows with no phenotype
    phenotypes = phenotypes.dropna(subset=["phenotype"])

    pheno_clf = []
    unknown_is_resistant = snakemake.params.unknown_is_resistant
    failed_is_resistant = snakemake.params.failed_is_resistant
    classifier = Classifier(
        unknown_is_resistant=unknown_is_resistant,
        failed_is_resistant=failed_is_resistant,
    )

    for _, row in df.iterrows():
        drug = row["drug"]
        sample = row["sample"]
        w = row["w"]
        k = row["k"]
        tech = row["technology"]

        ph = phenotypes["phenotype"].get((sample, drug))
        if ph is None:
            continue
        truth = Prediction(ph)

        pred = Prediction(row["prediction"])
        clf = classifier.from_predictions(truth, pred)

        pheno_clf.append((sample, drug, str(clf), tech, w, k))

    clf_df = pd.DataFrame(
        pheno_clf,
        columns=[
            "sample",
            "drug",
            "classification",
            "technology",
            "w",
            "k",
        ],
    )

    wks = set(clf_df.loc[:, ["w", "k"]].itertuples(index=False, name=None))
    techs = set(clf_df["technology"])
    cms = []

    for t, (w, k) in product(techs, wks):
        s = clf_df.query("technology==@t and w==@w and k==@k").value_counts(
            subset=["classification"]
        )
        cms.append((t, w, k, *ConfusionMatrix.from_series(s).ravel()))

    data = pd.DataFrame(cms, columns=["technology", "w", "k", "TN", "FP", "FN", "TP"])

    # set aesthetics
    plt.style.use(snakemake.params.style)
    plt.rcParams["figure.figsize"] = snakemake.params.figsize
    plt.rcParams["figure.dpi"] = snakemake.params.dpi

    x = "k"
    hue = "w"
    fs = snakemake.params.fontsize
    linestyles = list(lines.lineStyles.keys())[: len(set(data[hue]))]
    pointmarkers = ["o", "^"]
    fig, axes = plt.subplots(ncols=2, sharex=True, sharey=True, squeeze=True)

    for i, tech in enumerate(techs):
        d = data.query("technology==@tech")
        ax = axes[i]
        ax.set_title(f"{tech.capitalize()}")

        handles = [Line2D([0], [0], color="none")]
        labels = [""]

        for j, y in enumerate(["FN", "FP"]):
            ls = linestyles[j]
            m = pointmarkers[j]

            sns.pointplot(
                data=d,
                x=x,
                y=y,
                hue=hue,
                ax=ax,
                dodge=True,
                ci=None,
                linestyles=ls,
                markers=m,
            )

            line = Line2D([], [], label=y, color="black", ls=ls, linewidth=2)
            point = plt.scatter([], [], s=100, facecolors="black", marker=m)
            handles.append((line, point))
            labels.append(y)

        if i == 0:
            ax.set_ylabel("FN/FP count", fontdict=dict(fontsize=fs))
        else:
            ax.set_ylabel("")

        leghandles, leglabels = ax.get_legend_handles_labels()
        leghandles = leghandles[: len(set(df[hue]))]
        leglabels = leglabels[: len(set(df[hue]))]
        leghandles.extend(handles)
        leglabels.extend(labels)

        ax.legend(
            leghandles,
            leglabels,
            loc="best",
            handlelength=3,
            fontsize=fs,
            title=hue,
            title_fontsize=fs,
        )
        ax.tick_params(axis="both", labelsize=fs)
        ax.xaxis.label.set_size(fs)

    plt.tight_layout()

    fig.savefig(snakemake.output.plot)

    data["technology"] = data["technology"].apply(str.capitalize)
    data.sort_values(by=["w", "k"]).to_csv(snakemake.output.table, index=False)


main()
