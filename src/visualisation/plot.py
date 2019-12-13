import os
import math
from functools import partial

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import confusion_matrix

from src.bio.util import subdirs


FILES = [
    ("metrics.csv", "epoch"),
    ("roc.csv", "index"),
    ("precision_recall.csv", "index"),
    ("average_precision.csv", "index"),
    ("auc.csv", "index"),
    ("predictions.csv", None),
]


def rgb(r, g, b):
    return r / 255.0, g / 255.0, b / 255.0


paletteG = partial(sns.light_palette, rgb(0, 61, 100), reverse=True, input="rgb")
gradientPalette = paletteG()
cmap = paletteG(as_cmap=True)
cmapI = paletteG(as_cmap=True, reverse=False)

sns.set_style("darkgrid")
# plt.rcParams.update({'font.size': 11})
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Source Sans Pro"  # ['Fira Sans', 'Source Sans Pro']
font = {"weight": "normal"}  # ,'size'   : 22}

palette = sns.color_palette(
    [
        rgb(0, 61, 100),  # UA Blue
        rgb(126, 0, 47),  # UA Red
        rgb(255, 140, 0),  # Orange
        rgb(153, 50, 204),  # Purple
        rgb(60, 204, 50),  # Green
        rgb(50, 201, 204),  # Cyan
        rgb(244, 66, 143),  # Pink
    ]
)

sns.set_palette(palette)
sns.set_style("whitegrid")


def getOutputPath(dir, title, extension=".pdf"):
    return os.path.join(dir, title + extension)


def consolidateAll(directory, force=False):
    # directory
    #   |- iteration 0
    #       |- metrics.csv
    #       |- roc.csv
    #       |- precision_recall.csv
    #       |- average_precision.csv
    #   |- iteration 1
    #       |- ...
    for file, col in FILES:
        outputPath = os.path.join(directory, file)

        if not force and os.path.exists(outputPath):
            continue

        consolidate(directory, file, col)


def consolidate(directory, file, col):
    dfs = list()
    for subdir in filter(
        lambda x: os.path.basename(x).startswith("iteration"), subdirs(directory)
    ):
        p = os.path.join(subdir, file)
        if not os.path.exists(p):
            print(f"{file} not found in one of the iterations, skipping.")
            return
        df = pd.read_csv(p)

        # Set roc end values to something that makes sense.
        if file == "roc.csv":
            df.tpr[0], df.tpr[-1] = 0.0, 1.0

        # Set p/r values to something that makes sense.
        if file == "precision_recall.csv":
            df.precision[0], df.recall[0] = 1.0, 0.0

        dfs.append(df)

    print(file)
    outputPath = os.path.join(directory, file)

    df_concat = pd.concat(dfs)
    if col is None:
        df_concat["type"] = os.path.basename(os.path.dirname(directory))
        df_concat.to_csv(outputPath, index=False)
        return
    elif col == "index":
        df_concat = df_concat.groupby(df_concat.index)
    elif col:
        df_concat = df_concat.groupby(df_concat[col], as_index=False)

    df_means = df_concat.mean()
    df_std = df_concat.std(ddof=0).add_prefix("std_")

    result = pd.concat([df_means, df_std], axis=1, sort=False)
    result["type"] = os.path.basename(os.path.dirname(directory))
    outputPath = os.path.join(directory, file)
    result.to_csv(outputPath, index=False)


def concatenateAll(directory, force=False):
    for file, col in FILES:
        # Can't concatenate unagregated csv's
        if col is None:
            continue

        outputPath = os.path.join(directory, file)

        if not force and os.path.exists(outputPath):
            continue

        concatenate(directory, file)


def concatenate(directory, file):
    dfs = list()
    for subdir in subdirs(directory):
        if os.path.basename(subdir).startswith("_"):
            print("Found directory:", subdir, "which starts with underscore, skipping")
            continue

        p = os.path.join(subdir, file)
        if not os.path.exists(p):
            print(f"{file} not found in one of the experiments, skipping.")
            return
        df = pd.read_csv(p)
        df["type"] = os.path.basename(subdir)
        dfs.append(df)

    df_concat = pd.concat(dfs)
    outputPath = os.path.join(directory, file)
    df_concat.to_csv(outputPath, index=False)


def plotMetrics(directory):
    metricsPath = os.path.join(directory, "metrics.csv")
    metrics = pd.read_csv(metricsPath)
    for metric in [
        "AUC",
        "acc",
        "balanced_accuracy",
        "loss",
        "mean_pred",
        "precision",
        "recall",
        "val_AUC",
        "val_acc",
        "val_balanced_accuracy",
        "val_loss",
        "val_mean_pred",
        "val_precision",
        "val_recall",
    ]:
        std_metric = "std_" + metric
        plt.figure()
        sns_plot = sns.lineplot(x="epoch", y=metric, ci=None, hue="type", data=metrics)

        labels = list()
        for tpe in metrics.type.unique():
            df = metrics[metrics.type == tpe]
            value = float(df.tail(1)[metric])
            labels.append("{} (final = {:.2f})".format(tpe, value))

        if metric in [
            "AUC",
            "acc",
            "balanced_accuracy",
            "precision",
            "recall",
            "val_AUC",
            "val_acc",
            "val_balanced_accuracy",
            "val_precision",
            "val_recall",
        ]:
            sns_plot.set_ylim(0.5, 1)
        elif metric in ["mean_pred", "val_mean_pred"]:
            sns_plot.set_ylim(0, 1)
        sns_plot.set_xlim(0,)
        # Force ticks to be ints
        sns_plot.xaxis.set_major_locator(MaxNLocator(integer=True))

        handles, _ = sns_plot.get_legend_handles_labels()
        sns_plot.legend(
            handles=handles[1:],
            labels=[l.capitalize() for l in labels],
            loc="upper left",
            title=None,
        )
        sns_plot.set_xlabel("Epoch")
        sns_plot.set_title(metric.capitalize())

        for tpe in metrics.type.unique():
            df = metrics[metrics.type == tpe]
            sns_plot.fill_between(
                df.epoch,
                df[metric] - df[std_metric],
                df[metric] + df[std_metric],
                alpha=0.5,
            )
        sns_plot.get_figure().savefig(
            getOutputPath(directory, metric), bbox_inches="tight"
        )


def plotRocPR(directory):

    # define facet figure
    fig = plt.figure(constrained_layout=True, dpi=200, figsize=(12, 6))
    gs = GridSpec(1, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # create roc plot
    plotRoc(directory, ax=ax1)

    # create precision-recall plot
    plotPrecisionRecall(directory, ax=ax2)

    # add labels
    import string

    for n, ax in enumerate(fig.axes):
        ax.text(
            -0.1,
            1.1,
            string.ascii_uppercase[n],
            transform=ax.transAxes,
            size=20,
            weight="bold",
        )

    fig.savefig(getOutputPath(directory, "roc-pr"), bbox_inches="tight")


def plotRoc(directory, ax=None):
    rocPath = os.path.join(directory, "roc.csv")
    roc = pd.read_csv(rocPath)

    aucPath = os.path.join(directory, "auc.csv")
    auc = pd.read_csv(aucPath)

    labels = list()
    for tpe in auc.type.unique():
        df = auc[auc.type == tpe]
        auc_mean = float(df.auc)
        auc_std = float(df.std_auc)
        labels.append("{} (AUC = {:.2f} ± {:.2f})".format(tpe, auc_mean, auc_std))

    save = False
    if ax is None:
        plt.figure(figsize=(6.4, 6.4))
        ax = plt.gca()
        save = True

    sns_plot = sns.lineplot(x="fpr", y="tpr", ci=None, data=roc, hue="type", ax=ax)

    sns_plot.set_title("ROC")
    sns_plot.legend(labels, title=None, loc="lower right")
    sns_plot.set_xlabel("False Positive Rate")
    sns_plot.set_ylabel("True Positive Rate")

    for tpe in roc.type.unique():
        df = roc[roc.type == tpe]
        sns_plot.fill_between(
            df.fpr, df.tpr - df.std_tpr, df.tpr + df.std_tpr, alpha=0.5
        )

    sns_plot.plot([0, 1], [0, 1], "k--")
    if save:
        sns_plot.get_figure().savefig(
            getOutputPath(directory, "roc"), bbox_inches="tight"
        )


def plotPrecisionRecall(directory, ax=None):
    precisionRecallPath = os.path.join(directory, "precision_recall.csv")
    precisionRecall = pd.read_csv(precisionRecallPath)

    # Interpolation messes these up if the highest predictions are negative samples.
    precisionRecall.at[0, "recall"] = 0
    precisionRecall.at[0, "precision"] = 1

    averagePrecisionPath = os.path.join(directory, "average_precision.csv")
    averagePrecision = pd.read_csv(averagePrecisionPath)

    labels = list()
    for tpe in averagePrecision.type.unique():
        df = averagePrecision[averagePrecision.type == tpe]
        prec_mean = float(df.average_precision)
        prec_std = float(df.std_average_precision)
        labels.append(
            "{} (Avg. Prec. = {:.2f} ± {:.2f})".format(tpe, prec_mean, prec_std)
        )

    save = False
    if ax is None:
        plt.figure(figsize=(6.4, 6.4))
        ax = plt.gca()
        save = True

    plt.figure(figsize=(6.4, 6.4))
    sns_plot = sns.lineplot(
        x="recall", y="precision", ci=None, data=precisionRecall, hue="type", ax=ax
    )

    sns_plot.set_title("Precision - Recall")
    sns_plot.legend(labels, title=None, loc="lower left")
    sns_plot.set_xlabel("Recall")
    sns_plot.set_ylabel("Precision")
    sns_plot.set_ylim(0, 1)

    for tpe in precisionRecall.type.unique():
        df = precisionRecall[precisionRecall.type == tpe]
        sns_plot.fill_between(
            df.recall,
            df.precision - df.std_precision,
            df.precision + df.std_precision,
            alpha=0.5,
        )

    if save:
        sns_plot.get_figure().savefig(
            getOutputPath(directory, "precision_recall"), bbox_inches="tight"
        )


def plotPredictions(directory):
    predictionsPath = os.path.join(directory, "predictions.csv")
    if not os.path.exists(predictionsPath):
        return

    predictions = pd.read_csv(predictionsPath)
    bins = np.linspace(0, 1, 41)
    plt.figure()
    sns_plot = sns.distplot(
        predictions.y_pred[predictions.y_true == 0], bins=bins, kde=False
    )
    sns_plot = sns.distplot(
        predictions.y_pred[predictions.y_true == 1], bins=bins, kde=False
    )
    sns_plot.set_xlim(0, 1)
    sns_plot.set_ylim(0, 18000)
    title = os.path.basename(os.path.normpath(directory))
    # title = "Predictions"
    sns_plot.set_title(title)
    sns_plot.set_xlabel("Predicted probability")
    sns_plot.legend(["Negative", "Positive"], title=None)
    sns_plot.get_figure().savefig(
        getOutputPath(directory, "predictions"), bbox_inches="tight"
    )


def plotConfusionMatrix(directory, ax=None):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    source: https://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html#sphx-glr-auto-examples-model-selection-plot-confusion-matrix-py
    """

    predictionsPath = os.path.join(directory, "predictions.csv")
    if not os.path.exists(predictionsPath):
        return

    predictions = pd.read_csv(predictionsPath)
    y_true = predictions.y_true
    y_pred = [1 if p > 0.5 else 0 for p in predictions.y_pred]

    classes = ["False", "True"]
    normalize = False

    save = False
    if ax is None:
        plt.figure()
        ax = plt.gca()
        save = True

    # Plot non-normalized confusion matrix
    np.set_printoptions(precision=2)

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)

    # Test samples are sampled 5 times. Normalize numbers back to size of validation set
    cm //= 5

    labels = [
        ["True Negatives", "False Positives"],
        ["False Negatives", "True Positives"],
    ]

    if normalize:
        cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]

    im = ax.imshow(cm, interpolation="nearest", cmap=cmapI)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(
        xticks=np.arange(0, cm.shape[1], 1),
        yticks=np.arange(0, cm.shape[0], 1),
        # ... and label them with the respective list entries
        xticklabels=classes,
        yticklabels=classes,
        title="Confusion Matrix",
        ylabel="Label",
        xlabel="Prediction",
    )

    # Minor ticks to plot grid
    ax.set_xticks(np.arange(-0.5, cm.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, cm.shape[0], 1), minor=True)

    ax.grid(False)
    ax.grid(which="minor")

    # Rotate the tick labels and set their alignment.
    # plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
    #          rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = ".2f" if normalize else "d"
    thresh = cm.max() / 2.0
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(
                j,
                i,
                f"{labels[i][j]}\n{format(cm[i, j], fmt)}",
                ha="center",
                va="center",
                color="white" if cm[i, j] > thresh else "black",
            )

    if save:
        ax.get_figure().savefig(
            getOutputPath(directory, "confusion_matrix"), bbox_inches="tight"
        )


def plotAll(directory):
    plotMetrics(directory)
    plotRoc(directory)
    plotPrecisionRecall(directory)
    plotRocPR(directory)
    plotPredictions(directory)
    plotConfusionMatrix(directory)


def plotCombined(directories):
    plotCombinedFunction(directories, plotRoc, "ROC")
    plotCombinedFunction(directories, plotPrecisionRecall, "Precision - Recall")


def plotCombinedFunction(directories, plotFunc, title):
    cols = 2
    rows = math.ceil(len(directories) / 2)
    width = 6.4 * cols
    height = 6.4 * rows
    f, axarr = plt.subplots(
        rows, cols, figsize=(width, height), sharey="none", sharex="none"
    )
    f.suptitle(title, y=0.91, fontsize=20)
    for index, directory in enumerate(directories):
        ax = axarr[index // cols][index % cols]
        plotFunc(directory, ax=ax)

        subtitle = os.path.basename(os.path.normpath(directory))
        ax.set_title(subtitle)

    f.savefig(getOutputPath("output", title), bbox_inches="tight")
