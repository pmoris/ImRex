from functools import partial
import math
import os
from pathlib import Path

# from textwrap import fill

import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from sklearn.metrics import (
    auc,
    average_precision_score,
    confusion_matrix,
    precision_recall_curve,
    roc_curve,
)

from src.bio.util import subdirs


FILES = [
    ("metrics.csv", "epoch"),
    ("roc.csv", "index"),
    ("auc.csv", "index"),
    ("precision_recall.csv", "index"),
    ("average_precision.csv", "index"),
    ("predictions.csv", None),
]


def rgb(r, g, b):
    return r / 255.0, g / 255.0, b / 255.0


palette_g = partial(sns.light_palette, rgb(0, 61, 100), reverse=True, input="rgb")
gradient_palette = palette_g()
cmap = palette_g(as_cmap=True)
cmap_i = palette_g(as_cmap=True, reverse=False)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Source Sans Pro"  # ['Fira Sans', 'Source Sans Pro']
plt.rcParams.update({"font.size": 16})
plt.rcParams["patch.edgecolor"] = "black"

sns.set_palette("Dark2")
sns.set_style("darkgrid")


def get_output_path(directory, title, extension=".pdf"):
    return os.path.join(directory, title + extension)


def derive_metrics_all(directory, force=False):
    """
    For each iteration, derive the roc and p/r values. Creates four new files, as shown below.

    directory
      |- iteration 0
          |- metrics.csv
          |- predictions.csv
          |+ roc.csv
          |+ auc.csv
          |+ precision_recall.csv
          |+ average_precision.csv
      |- iteration 1
          |- ...
    """
    for subdir in filter(
        lambda x: os.path.basename(x).startswith("iteration"), subdirs(directory)
    ):
        p = os.path.join(directory, os.path.basename(subdir))
        predictions_path = os.path.join(p, "predictions.csv")
        if not os.path.exists(predictions_path):
            print(f"Missing 'predictions.csv' in {p}...")
            continue
        if not os.path.getsize(predictions_path) > 0:
            print(f"{predictions_path} in {p} appears to be empty, skipping...")
            continue

        # read predictions from csv
        predictions = pd.read_csv(predictions_path, sep=",")
        y_pred, y_true = predictions.y_pred, predictions.y_true

        derive_roc(p, y_true, y_pred, force=force)
        derive_pr(p, y_true, y_pred, force=force)

        plot_predictions(subdir)


def derive_roc(subdir, y_true, y_pred, force=False):
    roc_path = os.path.join(subdir, "roc.csv")
    auc_path = os.path.join(subdir, "auc.csv")

    # Check if already processed
    if not force and os.path.exists(roc_path) and os.path.exists(auc_path):
        return

    # calculate fpr, tpr and thresholds for ROC curve
    fpr, tpr, _ = roc_curve(y_true, y_pred, drop_intermediate=True)

    # interpolate (to easily combine with other iterations)
    interval = np.linspace(0, 1, 201)
    tpr_i = np.interp(interval, fpr, tpr)

    # Set roc end values to something that makes sense (should mostly be the case, up to rounding errors).
    tpr_i[0], tpr_i[-1] = 0.0, 1.0

    # write to file
    df = pd.DataFrame({"fpr": interval, "tpr": tpr_i})
    df.to_csv(roc_path, index=False)

    # calculate auc from ROC curve. This is done here already (and written to file) so it can be processed by the same average/stddev calculations later on.
    auc_value = auc(fpr, tpr)

    # write to file
    df = pd.DataFrame({"auc": [auc_value]})
    df.to_csv(auc_path, index=False)


def derive_pr(subdir, y_true, y_pred, force=False):
    pr_path = os.path.join(subdir, "precision_recall.csv")
    apr_path = os.path.join(subdir, "average_precision.csv")

    # Check if already processed
    if not force and os.path.exists(pr_path) and os.path.exists(apr_path):
        return

    # calculate p/r values
    precision, recall, _ = precision_recall_curve(y_true, y_pred)

    # interpolate (to easily combine with other iterations)
    interval = np.linspace(0, 1, 201)
    precision_i = np.interp(interval, recall[::-1], precision[::-1])

    # Set p/r values to something that makes sense (should mostly be the case, up to rounding errors).
    precision_i[0], interval[0] = 1.0, 0.0

    # write to file
    df = pd.DataFrame({"recall": interval, "precision": precision_i})
    df.to_csv(pr_path, index=False)

    # calculate average precitions. This is done here already (and written to file) so it can be processed by the same average/stddev calculations later on.
    ap = average_precision_score(y_true, y_pred)

    # write to file
    df = pd.DataFrame({"average_precision": [ap]})
    df.to_csv(apr_path, index=False)


def consolidate_all(directory, force=False):
    # directory
    #   |- iteration 0
    #       |- metrics.csv
    #       |- roc.csv
    #       |- precision_recall.csv
    #       |- average_precision.csv
    #   |- iteration 1
    #       |- ...
    for file, col in FILES:
        output_path = os.path.join(directory, file)

        if not force and os.path.exists(output_path):
            continue

        consolidate(directory, file, col)

    # create auc per epitope csv for later box plot auroc comparisons
    consolidate_auc(directory)


def consolidate_auc(directory):
    dfs = list()
    for subdir in filter(
        lambda x: os.path.basename(x).startswith("iteration"), subdirs(directory)
    ):
        p = os.path.join(subdir, "auc.csv")

        if not os.path.exists(p):
            print(f"auc.csv not found in {subdir}, skipping...")
            continue

        df = pd.read_csv(p)
        df["iteration"] = os.path.basename(subdir)
        dfs.append(df)

    if not dfs:
        print(f"No auc.csv found in any subdirectory of {directory}, skipping...")
        return

    output_path = os.path.join(directory, "auc_per_iteration.csv")
    df_concat = pd.concat(dfs)
    df_concat["type"] = os.path.basename(os.path.abspath(directory))
    df_concat.to_csv(output_path, index=False)


def consolidate(directory, file, col):
    dfs = list()
    for subdir in filter(
        lambda x: os.path.basename(x).startswith("iteration"), subdirs(directory)
    ):
        p = os.path.join(subdir, file)

        if not os.path.exists(p):
            print(f"{file} not found in {subdir}, skipping...")
            continue
        if not os.path.getsize(p) > 0:
            print(f"{file} in {subdir} appears to be empty, skipping...")
            continue

        df = pd.read_csv(p)

        # Moved to derive step. If still needed for legacy results -> uncomment
        # # Set roc end values to something that makes sense.
        # if file == "roc.csv":
        #     df.tpr[0], df.tpr[-1] = 0.0, 1.0
        #
        # # Set p/r values to something that makes sense.
        # if file == "precision_recall.csv":
        #     df.precision[0], df.recall[0] = 1.0, 0.0

        dfs.append(df)

    # skip if no files were found
    if not dfs:
        print(f"No {file} found in any subdirectory of {directory}, skipping...")
        return

    output_path = os.path.join(directory, file)
    df_concat = pd.concat(dfs)

    if col is None:
        df_concat["type"] = os.path.basename(os.path.abspath(directory))
        df_concat.to_csv(output_path, index=False)
        return
    elif col == "index":
        df_concat = df_concat.groupby(df_concat.index)
    elif col:
        df_concat = df_concat.groupby(df_concat[col], as_index=False)

    df_means = df_concat.mean()
    df_std = df_concat.std(ddof=0).add_prefix("std_")

    result = pd.concat([df_means, df_std], axis=1, sort=False)

    result["type"] = os.path.basename(os.path.abspath(directory))
    result.to_csv(output_path, index=False)


def concatenate_all(directory, force=False):
    for file, col in FILES + [("auc_per_iteration.csv", "index")]:
        # Can't concatenate unagregated csv's
        if col is None:
            continue

        output_path = os.path.join(directory, file)

        if not force and os.path.exists(output_path):
            continue

        concatenate(directory, file)


def concatenate(directory, file):
    dfs = list()
    for subdir in subdirs(directory):
        if os.path.basename(subdir).startswith("_"):
            print(
                f"Found directory: {subdir} which starts with underscore, skipping..."
            )
            continue

        p = os.path.join(subdir, file)
        if not os.path.exists(p):
            print(f"{file} not found in one of the experiments, skipping...")
            return

        df = pd.read_csv(p)
        df["type"] = os.path.basename(subdir)
        dfs.append(df)

    # skip if no files were found
    if not dfs:
        print(f"No {file} found in any subdirectory of {directory}, skipping...")
        return

    df_concat = pd.concat(dfs)
    output_path = os.path.join(directory, file)
    df_concat.to_csv(output_path, index=False)


def plot_metrics(directory, y_lim_loss=None):
    metrics_path = os.path.join(directory, "metrics.csv")

    if not os.path.exists(metrics_path):
        print(f"{metrics_path} not found, skipping plots...")
        return
    if not os.path.getsize(metrics_path) > 0:
        print(f"{metrics_path} appears to be empty, skipping plots...")
        return

    metrics = pd.read_csv(metrics_path)

    for metric in [
        "loss",
        "acc",
        "accuracy",
        "balanced_accuracy",
        "balanced_acc",
        "AUC",
        "auc",
        "roc_auc",
        "pr_auc",
        # "fn",
        # "fp",
        # "tn",
        # "tp",
        "precision",
        "recall",
        "val_loss",
        "val_acc",
        "val_accuracy",
        "val_balanced_accuracy",
        "val_balanced_acc",
        "val_AUC",
        "val_auc",
        "val_roc_auc",
        "val_pr_auc",
        # "val_fn",
        # "val_fp",
        # "val_tn",
        # "val_tp",
        "val_precision",
        "val_recall",
    ]:
        # check if metric is present, otherwise skip
        if metric not in metrics.columns:
            continue

        std_metric = "std_" + metric
        plt.figure()

        labels = list()
        if "type" in metrics.columns:
            for tpe in metrics.type.unique():
                df = metrics[metrics.type == tpe]
                value = float(df.tail(1)[metric])
                value_std = float(df.tail(1)[std_metric])
                labels.append(
                    "{}\n(final = {:.4g} ± {:.4g} ".format(tpe, value, value_std)
                    + r"$s$)"
                )

            sns_plot = sns.lineplot(
                x="epoch", y=metric, ci=None, hue="type", data=metrics
            )

            for tpe in metrics.type.unique():
                df = metrics[metrics.type == tpe]
                sns_plot.fill_between(
                    df.epoch,
                    df[metric] - df[std_metric],
                    df[metric] + df[std_metric],
                    alpha=0.5,
                )

        else:
            value = float(metrics.tail(1)[metric])
            value_std = float(metrics.tail(1)[std_metric])
            labels.append(
                "{}\n(final = {:.4g} ± {:.4g} ".format(directory, value, value_std)
                + r"$s$)"
            )
            sns_plot = sns.lineplot(x="epoch", y=metric, ci=None, data=metrics)

        handles, _ = sns_plot.get_legend_handles_labels()
        sns_plot.legend(
            handles=handles[1:],
            labels=[l.capitalize() for l in labels],
            loc="upper left",
            bbox_to_anchor=(1, 1),
            title=None,
        )

        if metric in [
            "acc",
            "accuracy",
            "balanced_accuracy",
            "AUC",
            "auc",
            "roc_auc",
            "pr_auc",
            "precision",
            "recall",
            "val_acc",
            "val_accuracy",
            "val_balanced_accuracy",
            "val_balanced_acc",
            "val_AUC",
            "val_auc",
            "val_roc_auc",
            "val_pr_auc",
            "val_precision",
            "val_recall",
        ]:
            sns_plot.set_ylim(0.5, 1)
        elif metric in ["mean_pred", "val_mean_pred"]:
            sns_plot.set_ylim(0, 1)
        elif metric in ["loss", "val_loss"] and y_lim_loss:
            sns_plot.set_ylim(df[metric].min() * 0.9, y_lim_loss)
        sns_plot.set_xlim(0,)
        # Force ticks to be ints
        sns_plot.xaxis.set_major_locator(MaxNLocator(integer=True))

        sns_plot.set_xlabel("Epoch")
        sns_plot.set_title(metric.capitalize())

        sns_plot.get_figure().savefig(
            get_output_path(
                directory,
                metric + "-" + str(Path(directory).absolute().name).replace(" ", "-"),
            ),
            bbox_inches="tight",
        )


def plot_loss(directory, y_lim_loss=None):
    metrics_path = os.path.join(directory, "metrics.csv")
    if not os.path.exists(metrics_path):
        print(f"{metrics_path} not found, skipping plots...")
        return
    if not os.path.getsize(metrics_path) > 0:
        print(f"{metrics_path} appears to be empty, skipping plots...")
        return

    metrics = pd.read_csv(metrics_path)

    if "loss" not in metrics.columns or "val_loss" not in metrics.columns:
        print("No loss values found in metrics.csv, skipping...")
        return

    train_loss_df = metrics.loc[:, ["loss", "std_loss", "epoch", "type"]]
    train_loss_df["train_val"] = "train"
    train_loss_df = train_loss_df.rename(
        columns={"loss": "value", "std_loss": "std_value"}
    )

    val_loss_df = metrics.loc[:, ["val_loss", "std_val_loss", "epoch", "type"]]
    val_loss_df["train_val"] = "val"
    val_loss_df = val_loss_df.rename(
        columns={"val_loss": "value", "std_val_loss": "std_value"}
    )

    loss_df = pd.concat([train_loss_df, val_loss_df])
    loss_df["type_train_val"] = loss_df["type"] + "_" + loss_df["train_val"]

    plt.figure()

    labels = list()
    if "type" in loss_df.columns:
        for tpe in loss_df.type_train_val.unique():
            df = loss_df.loc[loss_df.type_train_val == tpe]
            value = float(df.tail(1)["value"])
            value_std = float(df.tail(1)["std_value"])
            labels.append(
                "{}\n(final = {:.4g} ± {:.4g} ".format(tpe, value, value_std) + r"$s$)"
            )

        sns_plot = sns.lineplot(
            x="epoch", y="value", ci=None, hue="type_train_val", data=loss_df
        )

        for tpe in loss_df.type_train_val.unique():
            df = loss_df.loc[loss_df.type_train_val == tpe]
            sns_plot.fill_between(
                df.epoch,
                df["value"] - df["std_value"],
                df["value"] + df["std_value"],
                alpha=0.5,
            )

    else:
        for i in loss_df.train_val.unique():
            value = float(df.loc[df.train_val == i].tail(1)["value"])
            value_std = float(df.tail(1)["std_value"])
            labels.append(
                "{}\n(final = {:.4g} ± {:.4g} ".format(directory, value, value_std)
                + r"$s$)"
            )
        sns_plot = sns.lineplot(
            x="epoch", y="value", ci=None, hue="type_train_val", data=loss_df
        )

    handles, _ = sns_plot.get_legend_handles_labels()
    sns_plot.legend(
        handles=handles[1:],
        labels=[l.capitalize() for l in labels],
        loc="upper left",
        bbox_to_anchor=(1, 1),
        title=None,
    )

    sns_plot.set_xlim(0,)
    if y_lim_loss:
        sns_plot.set_ylim(loss_df["value"].min() * 0.9, y_lim_loss)
    # Force ticks to be ints
    sns_plot.xaxis.set_major_locator(MaxNLocator(integer=True))

    sns_plot.set_xlabel("Epoch")
    sns_plot.set_title("Train and validation loss")

    sns_plot.get_figure().savefig(
        get_output_path(
            directory,
            "train_val_loss"
            + "-"
            + str(Path(directory).absolute().name).replace(" ", "-"),
        ),
        bbox_inches="tight",
    )


def plot_roc_pr(directory):

    # define facet figure
    fig = plt.figure(
        constrained_layout=True,
        dpi=200,
        figsize=(16, 8),  # (12, 6),  # (14,8) is better for overlapping legends
    )
    gs = GridSpec(1, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # create roc plot
    roc_success = plot_roc(directory, ax=ax1, legend=True)

    # create precision-recall plot
    pr_success = plot_precision_recall(directory, ax=ax2, legend=True)

    # skip if either roc or pr are not present
    if not roc_success or not pr_success:
        return

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

    # for ax in fig.axes:
    #     plt.setp(
    #         ax.get_legend().get_texts(), fontsize="10.5"
    #     )  # "13")  # for legend text
    #     # plt.setp(ax.get_legend().get_title(), fontsize="11")  # "13")
    #     # ax.axis("equal")

    fig.savefig(
        get_output_path(
            directory,
            "roc-pr" + "-" + str(Path(directory).absolute().name).replace(" ", "-"),
        ),
        bbox_inches="tight",
    )


def plot_roc(directory, ax=None, legend=True):
    roc_path = os.path.join(directory, "roc.csv")

    if not os.path.exists(roc_path):
        print(f"{roc_path} not found, skipping roc plot...")
        return
    if not os.path.getsize(roc_path) > 0:
        print(f"{roc_path} appears to be empty, skipping roc plot...")
        return

    roc = pd.read_csv(roc_path)

    auc_path = os.path.join(directory, "auc.csv")
    auc = pd.read_csv(auc_path)

    save = False
    if ax is None:
        plt.figure(figsize=(6.4, 6.4))
        ax = plt.gca()
        save = True

    labels = list()
    if "type" in auc.columns:
        for tpe in auc.type.unique():
            df = auc[auc.type == tpe]
            auc_mean = float(df.auc)
            auc_std = float(df.std_auc)
            labels.append(
                f"{tpe}\n"
                + r"($\overline{AUC}$"
                + " = {:.4g} ± {:.4g} ".format(auc_mean, auc_std)
                + r"$s$)"
            )
            # labels = [fill(l, 60) for l in labels]
        sns_plot = sns.lineplot(
            x="fpr",
            y="tpr",
            ci=None,
            data=roc,
            hue="type",
            ax=ax,
            # palette=get_palette(roc, "type"),
        )

    else:
        auc_mean = float(auc.auc)
        labels.append(
            "{}\n(AUC = {:.4g})".format(os.path.basename(directory), auc_mean)
        )
        # labels = [fill(l, 60) for l in labels]
        sns_plot = sns.lineplot(x="fpr", y="tpr", ci=None, data=roc, ax=ax)

    sns_plot.set_title("ROC")

    if legend:
        sns_plot.legend(
            labels, title=None, loc="upper center", bbox_to_anchor=(0.5, -0.15)
        )  # loc="lower right")

    sns_plot.set_xlabel("False Positive Rate")
    sns_plot.set_ylabel("True Positive Rate")
    # sns_plot.set_ylim(-0.05, 1.05)
    # sns_plot.set_xlim(-0.05, 1.05)

    if "type" in roc.columns:
        for tpe in roc.type.unique():
            df = roc[roc.type == tpe]
            sns_plot.fill_between(
                df.fpr, df.tpr - df.std_tpr, df.tpr + df.std_tpr, alpha=0.5
            )

    sns_plot.plot([0, 1], [0, 1], "k--")
    if save:
        sns_plot.get_figure().savefig(
            get_output_path(
                directory,
                "roc" + "-" + str(Path(directory).absolute().name).replace(" ", "-"),
            ),
            bbox_inches="tight",
        )
    return True


def plot_precision_recall(directory, ax=None, legend=True):
    precision_recall_path = os.path.join(directory, "precision_recall.csv")

    if not os.path.exists(precision_recall_path):
        print(f"{precision_recall_path} not found, skipping precision recall plot...")
        return
    if not os.path.getsize(precision_recall_path) > 0:
        print(
            f"{precision_recall_path} appears to be empty, skipping precision recall plot..."
        )
        return

    precision_recall = pd.read_csv(precision_recall_path)

    # Interpolation messes these up if the highest predictions are negative samples.
    precision_recall.at[0, "recall"] = 0
    precision_recall.at[0, "precision"] = 1

    average_precision_path = os.path.join(directory, "average_precision.csv")
    average_precision = pd.read_csv(average_precision_path)

    save = False
    if ax is None:
        plt.figure(figsize=(6.4, 6.4))
        ax = plt.gca()
        save = True

    plt.figure(figsize=(6.4, 6.4))

    labels = list()
    if "type" in average_precision.columns:
        for tpe in average_precision.type.unique():
            df = average_precision[average_precision.type == tpe]
            prec_mean = float(df.average_precision)
            prec_std = float(df.std_average_precision)
            labels.append(
                "{}\n(Avg. Prec. = {:.4g} ± {:.4g} ".format(tpe, prec_mean, prec_std)
                + r"$s$)"
            )
            # labels = [fill(l, 60) for l in labels]
        sns_plot = sns.lineplot(
            x="recall", y="precision", ci=None, data=precision_recall, hue="type", ax=ax
        )
    else:
        prec_mean = float(average_precision.average_precision)
        labels.append(
            "{}\n(Avg. Prec. = {:.4g})".format(os.path.basename(directory), prec_mean)
        )
        # labels = [fill(l, 60) for l in labels]
        sns_plot = sns.lineplot(
            x="recall", y="precision", ci=None, data=precision_recall, ax=ax
        )

    sns_plot.set_title("Precision - Recall")

    if legend:
        sns_plot.legend(
            labels, title=None, loc="upper center", bbox_to_anchor=(0.5, -0.15)
        )

    # sns_plot.legend(labels, title=None, loc="lower left")
    sns_plot.set_xlabel("Recall")
    sns_plot.set_ylabel("Precision")
    # sns_plot.set_ylim(-0.05, 1.05)
    # sns_plot.set_xlim(-0.05, 1.05)

    if "type" in precision_recall.columns:
        for tpe in precision_recall.type.unique():
            df = precision_recall[precision_recall.type == tpe]
            sns_plot.fill_between(
                df.recall,
                df.precision - df.std_precision,
                df.precision + df.std_precision,
                alpha=0.5,
            )

    if save:
        sns_plot.get_figure().savefig(
            get_output_path(
                directory,
                "precision_recall"
                + "-"
                + str(Path(directory).absolute().name).replace(" ", "-"),
            ),
            bbox_inches="tight",
        )
    return True


def plot_predictions(directory):
    predictions_path = os.path.join(directory, "predictions.csv")

    if not os.path.exists(predictions_path):
        print(f"{predictions_path} not found, skipping predictions plot...")
        return
    if not os.path.getsize(predictions_path) > 0:
        print(f"{predictions_path} appears to be empty, skipping predictions plot...")
        return

    predictions = pd.read_csv(predictions_path)
    bins = np.linspace(0, 1, 41)
    plt.figure()
    sns_plot = sns.distplot(
        predictions.y_pred[predictions.y_true == 0], bins=bins, kde=False
    )
    sns_plot = sns.distplot(
        predictions.y_pred[predictions.y_true == 1], bins=bins, kde=False
    )
    sns_plot.set_xlim(0, 1)
    title = os.path.basename(os.path.normpath(os.path.abspath(directory)))
    # title = "Predictions"
    sns_plot.set_title(title)
    sns_plot.set_xlabel("Predicted probability")
    sns_plot.legend(["Negative", "Positive"], title=None)
    sns_plot.get_figure().savefig(
        get_output_path(
            directory,
            "predictions"
            + "-"
            + str(Path(directory).absolute().name).replace(" ", "-"),
        ),
        bbox_inches="tight",
    )


def plot_confusion_matrix(directory, ax=None):
    """ Print and plot the confusion matrix.

    Normalization can be applied by setting `normalize=True`.
    source: https://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html#sphx-glr-auto-examples-model-selection-plot-confusion-matrix-py
    """
    predictions_path = os.path.join(directory, "predictions.csv")

    if not os.path.exists(predictions_path):
        print(f"{predictions_path} not found, skipping confusion matrix plot...")
        return
    if not os.path.getsize(predictions_path) > 0:
        print(
            f"{predictions_path} appears to be empty, skipping confusion matrix plot..."
        )
        return

    predictions = pd.read_csv(predictions_path)
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

    im = ax.imshow(cm, interpolation="nearest", cmap=cmap_i)
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
    fmt = ".4g" if normalize else "d"
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
            get_output_path(
                directory,
                "confusion_matrix"
                + "-"
                + str(Path(directory).absolute().name).replace(" ", "-"),
            ),
            bbox_inches="tight",
        )


def plot_roc_boxplot(directory):
    fig, ax = plt.subplots(
        constrained_layout=False, dpi=200, figsize=(12, 6)
    )  # figsize=(12, 6), figsize=(14, 8) for large comparisons
    df = pd.read_csv(os.path.join(directory, "auc_per_iteration.csv"))

    # labels = list()
    for tpe in df.type.unique():
        df_label = df[df.type == tpe]
        auc_mean = df_label.auc.mean()
        auc_std = df_label.auc.std()
        model_name = (
            f"{tpe}\n"
            + r"($\overline{AUROC}$"
            + " = {:.4g} ± {:.4g} ".format(auc_mean, auc_std)
            + r"$s$)"
        )
        # labels.append(model_name)
        # labels = [fill(l, 50) for l in labels]

        df.loc[df.type == tpe, "type-mean-std"] = model_name

    sns_plot = sns.boxplot(
        x="type-mean-std",
        y="auc",
        # data=df,  # .sort_values(by=["auc"], ascending=False),
        # order=df.sort_values(by=["auc"], ascending=False),
        data=df.sort_values(by=["type-mean-std"]),
        order=sorted(df["type-mean-std"].unique()),
        palette=get_palette(df, "type-mean-std"),
        hue="type-mean-std",
    )
    # hue="type-mean-std" for legend, optionally use custom labels

    plt.setp(
        ax.get_xticklabels(), rotation=45, va="top", rotation_mode="default",
    )

    plt.xlabel(None)
    plt.ylabel("AUROC")

    # for patch in ax.artists:
    #     r, g, b, a = patch.get_facecolor()
    #     patch.set_facecolor((r, g, b, 0.7))
    # plt/sns_plot.legend() overrides custom legend
    # for lh in plt.legend().legendHandles:
    #     lh.set_alpha(0.7)

    # add legend and remove x labels
    # must be called after setting alpha or it will override location again
    ## ax.legend()

    # legend below figure, requires hue to be set
    sns_plot.legend(title=None, loc="upper center", bbox_to_anchor=(0.5, -0.15))
    sns_plot.set(xticklabels=[])

    # add wmu test

    if df["type"].nunique() == 2:
        # sort the values
        df = df.sort_values(["iteration", "type"]).reset_index(drop=True)
        # and calculate the difference in auroc between the model types per iteration
        df["diff"] = df.groupby("iteration")["auc"].transform(lambda x: x.diff())

        # # wilcoxon signed rank test
        # # compute wilcoxon test and add p-value to legend
        # p = scipy.stats.wilcoxon(df["diff"].dropna(), correction=True, mode="exact")[1]
        # ax.legend(title=r"Wilcoxon signed-rank test $P-value =$ " + str(round(p, 4)),)

        # # sign test
        # x = (df["diff"].dropna() > 0).sum()
        # n = len(df["diff"].dropna())
        # p = scipy.stats.binom_test(x, n, 0.5, alternative="two-sided")
        # ax.legend(title=f"Sign test {p}")

        # MWU test
        auroc_lists = [df.loc[df["type"] == i, "auc"] for i in df["type"].unique()]
        p = scipy.stats.mannwhitneyu(auroc_lists[0], auroc_lists[1])[1]

        sns_plot.legend(
            title="MWU test P-value = " + "{:0.4g}".format(p),  # {round(p,4)}",
            loc="upper center",
            bbox_to_anchor=(0.5, -0.15),
        )

    sns_plot.get_figure().savefig(
        get_output_path(
            directory,
            "roc_boxplot"
            + "-"
            + str(Path(directory).absolute().name).replace(" ", "-"),
        ),
        bbox_inches="tight",
    )


def plot_all(directory, y_lim_loss=None):
    plot_metrics(directory, y_lim_loss=y_lim_loss)
    plt.close("all")
    plot_loss(directory, y_lim_loss=y_lim_loss)
    plt.close("all")
    plot_roc(directory)
    plt.close("all")
    plot_precision_recall(directory)
    plt.close("all")
    plot_roc_pr(directory)
    plt.close("all")
    plot_predictions(directory)
    plt.close("all")
    plot_confusion_matrix(directory)
    plt.close("all")
    plot_roc_boxplot(directory)
    plt.close("all")


def plot_combined(directories):
    plot_combined_function(directories, plot_roc, "ROC")
    plot_combined_function(directories, plot_precision_recall, "Precision - Recall")


def plot_combined_function(directories, plot_func, title):
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
        plot_func(directory, ax=ax)

        subtitle = os.path.basename(os.path.normpath(os.path.abspath(directory)))
        ax.set_title(subtitle)

    f.savefig(
        get_output_path(
            "output",
            title + "-" + str(Path(directory).absolute().name).replace(" ", "-"),
        ),
        bbox_inches="tight",
    )


def roc_per_epitope(
    eval_df: pd.DataFrame,
    output_path: str,
    min_iterations: int = 5,
    comparison: bool = False,
    grouped: bool = False,
    decoy: bool = False,
):

    # avoid settingwithcopy warnings
    eval_df = eval_df.copy()

    plt.rcParams["patch.edgecolor"] = "black"

    # when comparing different model types
    if comparison:
        # only include epitopes that occurred in at least m iterations (= folds within a model)
        eval_df = eval_df[
            eval_df.groupby(["epitope", "type"])["epitope"].transform("count")
            >= min_iterations
        ]
        hue = "type"

        # omit all epitopes for which the highest auroc is lower than 0.5
        eval_df = eval_df[
            eval_df.groupby("epitope")["roc_auc"].transform(lambda x: max(x) > 0.5)
        ]
        # same as: eval_df.groupby("epitope")["roc_auc"].transform("max") >= 0.5

        # only include epitopes that have m iterations in each of the different models
        if not decoy:
            n_types = eval_df.type.nunique()
            eval_df = eval_df[
                eval_df.groupby("epitope").type.transform(
                    lambda x: x.nunique() == n_types
                )
            ]

        # for comparison with decoy models, only include epitopes that have m iterations in all decoy models
        # or m iterations in all normal models, i.e. "per epitope-set model type groups"
        else:
            n_types = eval_df.type.unique()
            n_types_decoy = len([i for i in n_types if "decoy" in i])
            n_types_normal = len(n_types) - n_types_decoy

            eval_df_decoy = eval_df[eval_df["type"].str.contains("decoy")]
            eval_df_normal = eval_df[~eval_df["type"].str.contains("decoy")]

            eval_df_decoy = eval_df_decoy[
                eval_df_decoy.groupby("epitope").type.transform(
                    lambda x: x.nunique() == n_types_decoy
                )
            ]

            eval_df_normal = eval_df_normal[
                eval_df_normal.groupby("epitope").type.transform(
                    lambda x: x.nunique() == n_types_normal
                )
            ]

            eval_df = pd.concat([eval_df_decoy, eval_df_normal])

        # # temp: only use top k roc_aucs. If a given epitope is selected because it's in the top k, also select it in the other model, regardless of its auroc
        # eval_df["max_roc_auc"] = eval_df.groupby(["epitope"])["roc_auc"].transform(
        #     "max"
        # )
        # top = 40
        # eval_df = (
        #     eval_df.sort_values(["max_roc_auc", "roc_auc"], ascending=False)
        #     .drop("max_roc_auc", axis=1)
        #     .reset_index(drop=True)
        #     .iloc[: 2 * top]
        # )

    # when evaluating a single model, only include epitopes that occurred in at least m iterations (= folds within the model)
    else:
        eval_df = eval_df[
            eval_df.groupby("epitope")["epitope"].transform("count") >= min_iterations
        ]
        hue = None

    fig, ax = plt.subplots(
        constrained_layout=True, dpi=200, figsize=(18, 6)  # (14, 6)
    )  # (16, 8))

    plotter = sns.barplot if grouped else sns.boxplot

    # if there is only 1 model, i.e. not a comparison plot
    # sort by roc only and use a single colour
    if not comparison:  # eval_df.type.nunique() == 1:
        order = (
            eval_df.groupby(["epitope"])["roc_auc"]
            .mean()
            .sort_values(ascending=False)
            .index
        )
        # sort dataframe according to order to simplify next steps
        # sorting must also happen on type in order to make colours match with legend
        eval_df["epitope"] = pd.Categorical(eval_df["epitope"], order)
        eval_df = eval_df.sort_values(["epitope"]).reset_index(drop=True)

        # if distance to training epitopes is provided
        # distance color bar is not used for kfold box plots because different iterations can have different distances
        if "min_dist" in eval_df.columns and grouped:

            dist = "min_dist"

            # # NOTE: colors are simply passed in their direct order to seaborn, i.e. the order argument is not utilised.
            # # consequently the colors should be sorted

            colors = mpl.cm.ScalarMappable(
                cmap=sns.light_palette(
                    sns.color_palette("Dark2")[0],
                    as_cmap=True,
                    # n_colors=eval_df[dist].nunique(),
                )
            ).to_rgba(eval_df[dist])

            cmap = sns.light_palette(
                sns.color_palette("Dark2")[0],
                as_cmap=True,
                # n_colors=eval_df[dist].nunique(),
            )
            norm = mpl.colors.BoundaryNorm(
                np.arange(eval_df[dist].min() - 0.5, eval_df[dist].max() + 1.5), cmap.N,
            )
            sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])

            cbar = plt.colorbar(
                sm, ticks=np.arange(eval_df[dist].min(), eval_df[dist].max() + 1)
            )
            cbar.set_label(
                "Minimum edit distance to training epitopes", rotation=270, labelpad=25,
            )

            g = plotter(
                x="epitope",
                y="roc_auc",
                # hue=hue,
                data=eval_df,
                # color=colour_palette,
                # palette=mpl.cm.ScalarMappable(cmap="magma").to_rgba(eval_df["mean_dist"]),
                palette=colors,
                order=order,
                # alpha=0.7,
            )

            # optionally overlap individual data points
            # sns.swarmplot(x="epitope", y="roc_auc", data=eval_df, color=".25")

            # create colours for distance in box plot manually
            # if not grouped and "min_dist" in eval_df.columns:
            #     eval_df = eval_df.sort_values("epitope")
            #     eval_df["color"] = (
            #         mpl.cm.ScalarMappable(cmap=cmap).to_rgba(eval_df[dist]).tolist()
            #     )

            #     for box, color in zip(
            #         g.artists, eval_df.drop_duplicates("epitope")["color"]
            #     ):
            #         box.set_facecolor(color)

        # if no distance is provided, or if non-grouped models are used:
        else:
            g = plotter(
                x="epitope",
                y="roc_auc",
                # hue=hue,
                data=eval_df,
                # color=colour_palette,
                # palette=mpl.cm.ScalarMappable(cmap="magma").to_rgba(eval_df["mean_dist"]),
                color=sns.color_palette("Dark2")[0],
                order=order,
                # alpha=0.7,
            )

    # comparison of multiple models
    else:
        # data = eval_df.sort_values(by="n", ascending=False)
        # data = eval_df.sort_values(by=["type", "roc_auc"], ascending=False)
        # order = eval_df.sort_values(by=["roc_auc"], ascending=False).index

        # sort by highest mean auroc across types
        order = (
            eval_df.groupby(["epitope", "type"])["roc_auc"]
            .mean()
            .sort_values(ascending=False)
            .groupby(level=0)
            .head(1)
            .index.get_level_values(0)
        )
        # sort by max/mean auroc
        # order = eval_df.groupby(["epitope"])["roc_auc"].max()

        # sort dataframe according to order to simplify next steps
        # sorting must also happen on type in order to make colours match with legend
        eval_df["epitope"] = pd.Categorical(eval_df["epitope"], order)
        eval_df = eval_df.sort_values(["epitope", "type"]).reset_index(drop=True)

        # single hue colors - in case no distance is provided
        colour_palette = get_palette(eval_df, "type")

        g = plotter(
            x="epitope",
            y="roc_auc",
            hue=hue,
            data=eval_df,
            palette=colour_palette,
            order=order,
            hue_order=eval_df.sort_values("type")["type"].unique()
            # edgecolor="black"
            # alpha=0.7
            # boxprops=dict(alpha=0.7),
        )

        # optionally use scatter plot instead of box plot
        # sns.scatterplot(
        #     x="mean_dist",
        #     y="roc_auc",
        #     hue="epitope",
        #     data=eval_df,
        #     style="type"
        #     # palette=colour_palette,
        # )

        # remove legend title (do not use for single model type plot, this will result in an empty rectangle)
        l = ax.legend()
        l.set_title("")

        # if len(colour_palette) == 7:
        #     # change alpha value of fill colours for custom palette, cannot be done through seaborn directly
        #     # see: https://github.com/mwaskom/seaborn/issues/979
        #     # add_alpha_to_legend(ax, l)
        #     for patch in ax.artists:
        #         r, g, b, a = patch.get_facecolor()
        #         patch.set_facecolor((r, g, b, 0.7))
        #     for lh in plt.legend().legendHandles:
        #         lh.set_alpha(0.7)

        # add wilcox test for non-decoy plots if there are exactly two model types, skip for decoy vs normal comparisons
        # for grouped => wilcox per epitope, for not-grouped (kfold) => wilcox per epitope mean
        if eval_df["type"].nunique() == 2 and not decoy:
            # calculate the difference in auroc between the model types, values must be sorted!
            if grouped:
                eval_df["diff"] = eval_df.groupby("epitope")["roc_auc"].transform(
                    lambda x: x.diff()
                )
                # compute wilcoxon test and add p-value to legend
                p = scipy.stats.wilcoxon(
                    eval_df["diff"].dropna(), correction=True  # , mode="exact"
                )[1]

                # print("Skew of differences")
                # print(eval_df["diff"].dropna())

                # print(
                #     [
                #         eval_df[eval_df["type"] == i].groupby("epitope")["roc_auc"].mean()
                #         for i in eval_df["type"].unique()
                #     ]
                # )
                print(scipy.stats.wilcoxon(eval_df["diff"].dropna(), correction=True))
                print(scipy.stats.normaltest(eval_df["diff"].dropna()))
                print(scipy.stats.shapiro(eval_df["diff"].dropna()))
                mean_roc_auc_lists = [
                    eval_df[eval_df["type"] == i].groupby("epitope")["roc_auc"].mean()
                    for i in eval_df["type"].unique()
                ]
                print(scipy.stats.ttest_rel(*mean_roc_auc_lists))

                # # sign test
                # x = (eval_df["diff"].dropna() > 0).sum()
                # n = len(eval_df["diff"].dropna())
                # p = scipy.stats.binom_test(x, n, 0.5, alternative="two-sided")  # [1]
                # ax.legend(title=f"Sign test {p}")

            # for non-grouped (i.e. k-fold), compare means of each epitope
            else:
                mean_roc_auc_lists = [
                    eval_df[eval_df["type"] == i].groupby("epitope")["roc_auc"].mean()
                    for i in eval_df["type"].unique()
                ]
                p = scipy.stats.wilcoxon(
                    *mean_roc_auc_lists, correction=True  # , mode="exact"
                )[1]
            ax.legend(
                title=r"Wilcoxon signed-rank test $P-value =$ " + "{:0.2g}".format(p),
                loc="upper right",
            )  # str(round(p, 4)),)

            # count number of epitopes for which one model is higher than the other, requires sort
            # eval_df.loc[eval_df["diff"] > 0,["epitope","type","roc_auc","diff"]].count()

        elif decoy and eval_df["type"].nunique() == 2:
            auroc_lists = [
                eval_df.loc[eval_df["type"] == i, "roc_auc"]
                for i in eval_df["type"].unique()
            ]
            p = scipy.stats.mannwhitneyu(auroc_lists[0], auroc_lists[1])[1]
            ax.legend(
                title=r"MWU test $P-value =$ " + "{:0.4g}".format(p)
            )  # str(round(p, 4)),)

        # color by training dist, must be positioned after legend creation to maintain main color instead of random first color gradient

        # set color bars if distance to training epitopes is provided
        # distance color bar is not used for kfold box plots because different iterations can have different distances
        # comparisons with decoy models are also skipped because it is difficult to keep track of which bar corresponds to which epitope/type
        dist = "min_dist"

        # if types differ in their distance range, set color bar + tick marks for every type, otherwise share them
        if (
            eval_df.groupby("type")[dist]
            .unique()
            .apply(lambda x: min(x))
            .unique()
            .shape[0]
            == 1
            and eval_df.groupby("type")[dist]
            .unique()
            .apply(lambda x: max(x))
            .unique()
            .shape[0]
            == 1
        ):
            shared_bar = True
            x_shift = 0.02
        else:
            shared_bar = False
            x_shift = 0.04

        if dist in eval_df.columns and grouped:
            # create color mappings per datapoint and color bar per type
            x_pos = 1.01
            n_types = eval_df["type"].nunique()
            n_epitopes = len(order)

            # must iterate through types in sorted order = same order as given to plot/legend
            for i, t in enumerate(eval_df.sort_values("type")["type"].unique()):

                # create color mapping
                eval_df.at[eval_df["type"] == t, "color"] = pd.Series(
                    mpl.cm.ScalarMappable(
                        cmap=sns.light_palette(
                            sns.color_palette("Dark2")[i], as_cmap=True
                        )
                    )
                    .to_rgba(eval_df.loc[eval_df["type"] == t, dist])
                    .tolist(),
                    index=eval_df.loc[eval_df["type"] == t, dist].index,
                )

                # for bar plots patches are created first for the first hue, then for the second, etc.
                # NOTE: The number of patches is equal for every type, i.e. even "absent" bars have an associated type.
                assert (
                    len(order) == len(ax.patches) // n_types
                ), "Mismatch between number of drawn bars and number of epitopes * directories."
                for j, epitope_bar in enumerate(order):

                    # sequentially retrieve the next patch, requires the plotting function to use sorted data
                    # e.g. both the legend, and the order of creating different bars, needs to use the same order of "type"
                    # otherwise this will not work for decoy models, because not every epitope will have the same types,
                    # so the order of the patches is difficult to predict (depends on order of occurrence in dataframe sorted by
                    # both colums simultaneously)
                    p = ax.patches[j + (i * n_epitopes)]
                    if (
                        epitope_bar
                        in eval_df.loc[eval_df["type"] == t, "epitope"].values
                    ):
                        # some epitopes might have multiple distances across iterations depending on train/test splits - test only required for non-grouped
                        # assert eval_df.loc[ (eval_df["type"] == t) & (eval_df["epitope"] == epitope_bar), dist].nunique() == 1, "Found multiple epitope distances for the same epitope across iterations."
                        color = eval_df.loc[
                            (eval_df["type"] == t)
                            & (eval_df["epitope"] == epitope_bar),
                            "color",
                        ].values[0]
                        p.set_color(color)
                        p.set_edgecolor("black")

                # OLD METHOD: fails when there is an inequal amount of rows in the dataframe per type (e.g. TRA evaluation)
                # yet the zip operation assumes that each patch has a corresponding row per type.
                # # color boxes / bars
                # if grouped:
                #     # patches are created first for the first hue, then for the second
                #     for p, color in zip(
                #         ax.patches[
                #             i
                #             * len(ax.patches)
                #             // n_types : (i + 1)
                #             * len(ax.patches)
                #             // n_types
                #         ],
                #         eval_df.loc[eval_df["type"] == t, "color"],
                #     ):
                #         p.set_color(color)
                #         p.set_edgecolor("black")
                # # else:
                #     # boxes are created alternating for hues per epitope
                #     for box, color in zip(
                #         g.artists[i::n],
                #         eval_df.drop_duplicates(["epitope", "type"]).loc[
                #             eval_df["type"] == t, "color"
                #         ],
                #     ):
                #         box.set_facecolor(color)

                # create color bar
                cmap = sns.light_palette(sns.color_palette("Dark2")[i], as_cmap=True)
                norm = mpl.colors.BoundaryNorm(
                    np.arange(
                        eval_df.loc[eval_df["type"] == t, dist].min() - 0.5,
                        eval_df.loc[eval_df["type"] == t, dist].max() + 1.5,
                    ),
                    cmap.N,
                )
                sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])

                # get bounding box information for the axes (since they're in a line, you only care about the top and bottom)
                bbox_ax = ax.get_position()

                # fig.add_axes() adds the colorbar axes
                # they're bounded by [x0, y0, x_width, y_width]
                cbar_im1a_ax = fig.add_axes(
                    [x_pos, bbox_ax.y0 + 0.07, 0.02, bbox_ax.y1 - bbox_ax.y0]
                )

                # only set tick marks on final color bar, unless the ranges differ
                if shared_bar:
                    cbar_im1a = plt.colorbar(sm, ticks=[], cax=cbar_im1a_ax)
                else:
                    cbar_im1a = plt.colorbar(
                        sm,
                        ticks=np.arange(
                            eval_df.loc[eval_df["type"] == t, dist].min(),
                            eval_df.loc[eval_df["type"] == t, dist].max() + 1,
                        ),
                        cax=cbar_im1a_ax,
                    )

                # increment color bar position for next one
                x_pos += x_shift

            # set tick marks for final color bar (overwrites the last empty one)
            if shared_bar:
                cbar_im1a = plt.colorbar(
                    sm,
                    ticks=np.arange(eval_df[dist].min(), eval_df[dist].max() + 1),
                    cax=cbar_im1a_ax,
                )

            # add label to final color bar
            cbar_im1a.set_label(
                "Minimum edit distance to training epitopes", rotation=270, labelpad=25
            )

    # rotate epitope labels
    plt.setp(
        ax.get_xticklabels(), rotation=90, va="top", rotation_mode="default",
    )
    # plt.setp(
    #     ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor",
    # )

    # set axis labels
    ax.set_ylabel("AUROC")
    ax.set_xlabel("Epitope")

    # set y-axis limit and format tick marks
    if grouped:
        ax.set(ylim=(0, 1))
        # plt.yticks(np.arange(0, 1.1, 0.1)) # does not work?
        ax.yaxis.set_ticks(np.arange(0, 1.1, 0.1))
    else:
        # ax.set_ylim(eval_df.roc_auc.min() * 0.9, 1)
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.arange(round(start, 1), 1 + 0.1, 0.1))
    # ax.yaxis.set_ticks(np.arange(start, end+.1, 0.1)) # DONT DO THIS, AXIS WILL BE WRONG
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%0.1f"))

    plt.savefig(output_path, bbox_inches="tight")


def roc_train_corr(eval_df, output_path):

    # # omit all epitopes for which the highest auroc is lower than 0.5
    # eval_df = eval_df[
    #     eval_df.groupby("epitope")["roc_auc"].transform(lambda x: max(x) > 0.5)
    # ]
    # # when evaluating a single model, only include epitopes that occurred in at least m iterations (= folds within the model)
    # eval_df = eval_df[eval_df.groupby("epitope")["epitope"].transform("count") >= 25]

    fig, ax = plt.subplots(constrained_layout=True, dpi=200, figsize=(16, 8))
    g = sns.jointplot(y="roc_auc", x="train_size", data=eval_df)
    # g.fig.subplots_adjust(top=0.93, wspace=0.3)
    # g.fig.suptitle("")
    # g.ax_joint.set_ylim((0,1))

    g.ax_joint.set_xlabel("Number of training observations (log-scale)")
    g.ax_joint.set_ylabel("AUROC")

    g.ax_joint.set_xscale("log")
    g.ax_marg_x.set_xscale("log")
    # g.ax_joint.set_xlim(
    #     [eval_df["train_size"].min() - 100, eval_df["train_size"].max() + 100]
    # )
    # g.ax_marg_x.set_xlim(
    #     [eval_df["train_size"].min() - 100, eval_df["train_size"].max() + 100]
    # )

    # r, p = scipy.stats.spearmanr(eval_df["roc_auc"], eval_df["train_size"])
    # annot_kws = {"prop": {"family": "monospace", "weight": "bold", "size": 15}}
    # phantom = g.ax_joint.plot([], [], linestyle="", alpha=0)
    # g.ax_joint.legend([phantom], ["r={:f}, p={:f}".format(r, p)], **annot_kws)
    # g.ax_joint.legend(["r={:f}, p={:f}".format(r, p)])

    g.fig.set_dpi(200)

    plt.savefig(output_path, bbox_inches="tight")


def roc_min_dist_box(eval_df, output_path):
    g = sns.boxplot(
        y="roc_auc", x="min_dist", data=eval_df, color=sns.color_palette("Dark2")[0],
    )
    # sns.swarmplot(y="roc_auc", x="min_dist", data=eval_df, color=".25")
    g.set_xlabel("Minimum edit distance")
    g.set_ylabel("AUROC")

    r, p = scipy.stats.spearmanr(eval_df["roc_auc"], eval_df["min_dist"])
    # annot_kws = {"prop": {"family": "monospace", "weight": "bold", "size": 15}}
    # phantom = g.plot([], [], linestyle="", alpha=0)
    # g.legend([phantom], ["r={:f}, p={:f}".format(r, p)], **annot_kws)
    g.set_title("r={:f}, p={:f}".format(r, p))

    plt.savefig(output_path, bbox_inches="tight")


def roc_avg_dist_corr(eval_df, output_path):
    g = sns.jointplot(y="roc_auc", x="mean_dist", data=eval_df)
    g.ax_joint.set_xlabel("Average edit distance")
    g.ax_joint.set_ylabel("AUROC")
    g.fig.set_dpi(200)
    r, p = scipy.stats.spearmanr(eval_df["roc_auc"], eval_df["mean_dist"])
    annot_kws = {"prop": {"family": "monospace", "weight": "bold", "size": 15}}
    (phantom,) = g.ax_joint.plot([], [], linestyle="", alpha=0)
    g.ax_joint.legend([phantom], ["r={:f}, p={:f}".format(r, p)], **annot_kws)
    # g.ax_joint.legend(spearman)
    plt.savefig(output_path, bbox_inches="tight")


def roc_min_dist_corr(eval_df, output_path):
    g = sns.jointplot(y="roc_auc", x="min_dist", data=eval_df)
    g.ax_joint.set_xlabel("Minimum edit distance")
    g.ax_joint.set_ylabel("AUROC")
    g.fig.set_dpi(200)
    plt.savefig(output_path, bbox_inches="tight")


def roc_median_dist_corr(eval_df, output_path):
    g = sns.jointplot(y="roc_auc", x="median_dist", data=eval_df)
    g.ax_joint.set_xlabel("Median edit distance")
    g.ax_joint.set_ylabel("AUROC")
    g.fig.set_dpi(200)
    plt.savefig(output_path, bbox_inches="tight")


def get_palette(df, value):
    unique_values = sorted(df[value].unique())

    # if len(unique_values) < 8:
    #     pal = custom_palette
    # else:
    #     pal = sns.color_palette("Set1", n_colors=len(unique_values))

    pal = sns.color_palette("Dark2", n_colors=len(unique_values))

    palette_dict = dict(zip(unique_values, pal))
    return palette_dict


# def add_alpha_to_legend(ax, l):
#     # change alpha value of fill colours, cannot be done through seaborn directly
#     # see: https://github.com/mwaskom/seaborn/issues/979
#     for patch in ax.artists:
#         r, g, b, a = patch.get_facecolor()
#         patch.set_facecolor((r, g, b, 0.7))
#     for lh in l.legendHandles:
#         lh.set_alpha(0.7)


def plot_cv_folds(df, output_path):
    fold_indices = [
        (df.index.difference(df[df["fold"] == i].index), df[df["fold"] == i].index)
        for i in df["fold"].unique()
    ]
    X = df.cdr3
    y = df.y
    groups = pd.factorize(df["antigen.epitope"])[0]
    n_splits = df.fold.nunique()
    lw = 30

    cmap_data = plt.cm.tab20c
    cmap_cv = plt.cm.coolwarm

    fig, ax = plt.subplots(figsize=(20, 8))

    # Generate the training/testing visualizations for each CV split
    for ii, (tr, tt) in enumerate(fold_indices):
        # Fill in indices with the training/test groups
        indices = np.array([np.nan] * len(X))
        indices[tt] = 1
        indices[tr] = 0

        # Visualize the results
        ax.scatter(
            range(len(indices)),
            [ii + 0.5] * len(indices),
            c=indices,
            marker="_",
            lw=lw,
            cmap=cmap_cv,
            vmin=-0.2,
            vmax=1.2,
        )

    # Plot the data classes and groups at the end
    ax.scatter(
        range(len(X)), [ii + 1.5] * len(X), c=y, marker="_", lw=lw, cmap=cmap_data
    )

    # create color map that cycles for groups
    color_dict = {
        i: cmap_data.colors[i % len(cmap_data.colors)] for i in np.unique(groups)
    }
    color_list = [color_dict[i] for i in groups]
    ax.scatter(range(len(X)), [ii + 2.5] * len(X), marker="_", lw=lw, color=color_list)

    # Formatting
    yticklabels = list(range(n_splits)) + ["class", "group"]
    ax.set(
        yticks=np.arange(n_splits + 2) + 0.5,
        yticklabels=yticklabels,
        xlabel="Sample index (sorted by epitope and class label)",
        ylabel="CV iteration",
        ylim=[n_splits + 2.2, -0.2],
        xlim=[0, len(X)],
    )
    ax.set_title("Cross-validation folds")  # , fontsize=15)

    ax.legend(
        [Patch(color=cmap_cv(0.8)), Patch(color=cmap_cv(0.02))],
        ["Testing set", "Training set"],
        loc=(1.02, 0.8),
    )
    # Make the legend fit
    plt.tight_layout()
    fig.subplots_adjust(right=0.7)

    plt.savefig(output_path, bbox_inches="tight")

    plt.close("all")
