import itertools
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns

import src.bacli as bacli
from src.bio.image import image_from_matrices, image_from_matrix, image_from_tensor
from src.bio.peptide_feature import (  # noqa: I101
    Charge,
    Hydrophilicity,
    Hydrophobicity,
    Mass,
    Surface,
    Flexibility,
    Transfer,
    # TCRexBasicity,
    # TCRexHelicity,
    AtchleyFactor1,
    AtchleyFactor2,
    AtchleyFactor3,
    AtchleyFactor4,
    AtchleyFactor5,
    parse_features,
    parse_operator,
    IsoelectricPoint,
)
from src.bio.util import subdirs
from src.config import PROJECT_ROOT
from src.definitions.amino_acid_properties import AMINO_ACIDS
from src.visualisation.plot import (
    cmap,
    concatenate_all,
    consolidate_all,
    derive_metrics_all,
    plot_all,
    # palette,
    plot_cv_folds,
    plot_combined,
    plot_roc_boxplot,
    roc_min_dist_corr,
    roc_median_dist_corr,
    roc_min_dist_box,
    roc_per_epitope,
    roc_train_corr,
)


# OUTPUT_DIR = "output/"
OUTPUT_DIR = PROJECT_ROOT / "reports/figures"
SCALE = 50


@bacli.command
def test(model_file: str):
    """ Debug function used to test commands out. """
    from tensorflow.keras.models import load_model

    model = load_model(model_file)
    print(model.layers[0].get_weights())


@bacli.command
def render(model_file: str):
    """ Render dot representation of neural network. """
    from tensorflow.keras.models import load_model
    from keras.utils.vis_utils import plot_model

    model = load_model(model_file, compile=False)
    plot_model(
        model,
        to_file=model_file + ".pdf",
        show_shapes=True,
        show_layer_names=True,
        rankdir="TB",
    )


@bacli.command
def activations(
    model_file: str,
    epitope: str = None,
    cdr3: str = None,
    layer: str = "conv2d_1",
    vsc: bool = False,
    features: str = "mass,hydrophob,charge,hydrophil",
    operator: str = "best",
):
    """ Display the activations of a model for a given input. """
    from tensorflow.keras.models import load_model
    from keract import get_activations, display_activations  # , display_heatmaps

    from src.processing.image_padding import ImagePadding
    from src.bio.feature_builder import CombinedPeptideFeatureBuilder

    if epitope is None or cdr3 is None:
        print("Supply epitope and cdr3")
        exit()

    model = load_model(model_file, compile=False)
    if vsc:
        model.summary(line_length=80)
        model = model.get_layer("sequential_1")
    model.summary(line_length=80)
    model.compile(loss="mse", optimizer="adam")

    features_list = parse_features(features)
    operator = parse_operator(operator)
    feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)
    image = feature_builder.generate_feature([cdr3, epitope])
    padded, _ = ImagePadding(None, 20, 13).transform((image, None), padValue=0)

    x = np.array((padded,))
    activations = get_activations(model, x, layer)
    display_activations(activations, cmap="gray", save=False)
    # display_heatmaps(activations, padded, save=False)


def make_input(epitope, cdr3):
    # TODO: implement
    #  generate input image for epitope and cdr3 sequence
    #  need to handle 'X' amino acid
    from src.bio.feature_builder import CombinedPeptideFeatureBuilder
    from src.processing.image_generator import ImageGenerator
    from src.processing.image_padding import ImagePadding
    from src.processing.data_stream import DataStream

    features = "hydrophob,isoelectric,mass,hydrophil,charge"
    operator = "absdiff"
    width = 20
    height = 11

    features_list = parse_features(features)
    operator = parse_operator(operator)
    feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)

    stream = DataStream([((cdr3, epitope), None)])

    pos_images = ImageGenerator(stream, feature_builder)
    pos_out = ImagePadding(pos_images, width, height, pad_value=0)

    sample = pos_out.get()[0]
    return sample


@bacli.command
def cam(model_file: str, epitope, cdr3):
    from tensorflow.keras.models import load_model
    from vis.visualization import visualize_cam

    model = load_model(model_file)

    # last layer is combination of previous one
    # originalLayer, title = get_image(epitope, cdr3, False, True, 'best')[-1]

    # TODO: keras-vis code is outdated, try snippet online ()

    x = make_input(epitope, cdr3)

    for modifier in [None, "guided", "relu"]:
        cam = visualize_cam(
            model,
            -1,
            None,
            x,
            penultimate_layer_idx=9,
            backprop_modifier=modifier,
            grad_modifier=None,
        )
        # print(cam)
        camLayer = image_from_tensor(cam, mode="RGB")

        if modifier is None:
            modifier = "vanilla"

        print(modifier)

        layers = [
            # (originalLayer, "input", "CMYK"),
            (camLayer, modifier, "jet")
        ]

        img2plot(layers, epitope, cdr3, "CAM.pdf")


def perturbe(variables, symbol, amount):
    """ Generate all possible combinations from list 'variables' where 'amount' places are set to 'symbol'. """
    variables = list(variables)

    for indices in itertools.combinations(range(len(variables)), amount):
        variation = variables.copy()
        for index in indices:
            variation[index] = symbol
        yield variation


@bacli.command
def predict_variations(model_file: str, epitope, cdr3, aa="G", perturbations: int = 1):
    from tensorflow.keras.models import load_model

    model = load_model(model_file)

    samples = list()

    # first sample is original
    samples.append((epitope, cdr3))

    print(f"Generating variations with up to {perturbations} changes")

    variables = list(epitope + cdr3)
    for amount in range(1, perturbations + 1):
        for variation in perturbe(variables, aa, amount):
            var_epitope = "".join(variation[: len(epitope)])
            var_cdr3 = "".join(variation[len(epitope) :])
            samples.append((var_epitope, var_cdr3))

    # TODO: Need to test what replacing with X does. It works, but I don't know which values are used.

    # generate input images for each sample
    x = np.array([make_input(*sample) for sample in samples])

    # get scores of model
    predictions = np.squeeze(model.predict_on_batch(x).numpy())

    # output results
    print("=============== All results ===============")
    for sample, prediction in zip(samples, predictions):
        print(sample[0], sample[1], prediction)

    print("=============== Statistics ===============")
    basePrediction = predictions[0]
    mostDifferent = sorted(
        zip(samples, predictions), key=lambda el: el[1] - basePrediction, reverse=True
    )

    print("Base prediction:", basePrediction)
    print("Most deviating (more positive):")
    for sample, prediction in mostDifferent[:5]:
        if prediction >= basePrediction:
            print("\t", sample[0], sample[1], prediction)

    print("Most deviating (more negative):")
    for sample, prediction in reversed(mostDifferent[-5:]):
        if prediction <= basePrediction:
            print("\t", sample[0], sample[1], prediction)

    return samples, predictions


@bacli.command
def summary(model_file: str):
    """ Print summary of neural network. """
    from tensorflow.keras.models import load_model

    model = load_model(model_file)
    model.summary(line_length=80)


def get_image(
    epitope: str = None,
    cdr3: str = None,
    tensor: bool = False,
    cmyk: bool = False,
    operator: str = "best",
):

    for pep in [epitope, cdr3]:
        print(pep)

    mode = "CMYK" if cmyk else "RGB"

    matrices = list()
    layers = list()

    # features = [Mass(), Hydrophobicity(), Charge()]
    features = [Mass(), Hydrophobicity(), Hydrophilicity()]
    if cmyk:
        features.append(IsoelectricPoint())

    operator = parse_operator(operator)

    for index, feature in enumerate(features):
        matrix = feature.image_matrix(cdr3, epitope, operator=operator)
        print(matrix)
        matrices.append(matrix)
        if tensor:
            pixels = feature.image_tensor(cdr3, epitope)
        else:
            pixels = feature.image_matrix(cdr3, epitope, operator=operator)

        img = image_from_matrix(pixels, mode=mode, index=index)
        layers.append((img, feature.name))

    img = image_from_matrices(*matrices, mode=mode)
    layers.append((img, "Combined"))
    return layers


@bacli.command
def peptide(
    epitope: str = None,
    cdr3: str = None,
    tensor: bool = False,
    cmyk: bool = False,
    operator: str = "absdiff",
):
    """ Render peptide images. """

    plt.rcParams.update({"font.size": 14})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Source Sans Pro"
    font = {"weight": "normal"}

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    layers = get_image(epitope, cdr3, tensor, cmyk, operator)

    img2plot(layers, epitope, cdr3, "amino-acid-map.pdf")


def img2plot(layers, epitope, cdr3, name):
    fig = plt.figure(figsize=(12, 4))
    axes = fig.subplots(1, len(layers), sharey=True)
    for ax in fig.get_axes():
        ax.set_yticks(range(len(cdr3)))
        ax.set_yticklabels(cdr3)
        ax.set_xticks(range(len(epitope)))
        ax.set_xticklabels(epitope)
        ax.set(ylabel="CDR3")
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.text(0.5, 0.08, "Epitope", ha="center", va="center")

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in fig.get_axes():
        ax.label_outer()

    for i, (layer, title) in enumerate(layers):
        sub = fig.get_axes()[i]
        sub.set_title(title)
        rgb_im = layer.convert("RGB")
        pix = np.array(rgb_im)
        sub.imshow(pix, origin="lower")
        sub.grid(False)
        for spine in plt.gca().spines.values():
            spine.set_visible(False)

    print(f"saving figure {OUTPUT_DIR}")
    plt.savefig(os.path.join(OUTPUT_DIR, name), bbox_inches="tight")


@bacli.command
def features():
    data = list()
    features = [
        Charge(),
        Hydrophobicity(),
        Hydrophilicity(),
        IsoelectricPoint(),
        Mass(),
        Surface(),
        Flexibility(),
        Transfer(),
        # TCRexBasicity(),
        # TCRexHelicity(),
        AtchleyFactor1(),
        AtchleyFactor2(),
        AtchleyFactor3(),
        AtchleyFactor4(),
        AtchleyFactor5(),
    ]
    # features = [Charge(), Hydrophobicity(), IsoelectricPoint(), Mass(), Hydrophilicity()]
    # features = [Charge(), Hydrophobicity(), IsoelectricPoint(), Mass(), Hydrophilicity(), TCRexBasicity(), TCRexHelicity(), TCRexHydrophobicity(), TCRexMutationStability()]
    for f in features:
        print(f.name)
        values = list()
        for aa in AMINO_ACIDS:
            values.append(f._calculate(aa))
            print(f"{aa};{f._calculate(aa)}")
        print()
        data.append(values)

    df = pd.DataFrame(zip(*data), columns=[f.name for f in features])

    def corrfunc(x, y, **kws):
        r, _ = stats.pearsonr(x, y)
        # s, _ = stats.spearmanr(x, y)
        ax = plt.gca()
        ax.annotate("r = {:.2f}".format(r), xy=(0.05, 0.95), xycoords=ax.transAxes)

    # sns.set(font_scale=1.1)
    sns.set_style("ticks", {"axes.grid": False})
    g = sns.PairGrid(df, height=2)

    g.map_lower(sns.regplot)  # , s=25)
    g.map_diag(plt.hist)
    g.map_upper(sns.kdeplot, cmap=cmap)
    # g.map_lower(sns.kdeplot, cmap="Blues_d")
    g.map_lower(corrfunc)
    path = os.path.join(OUTPUT_DIR, "features-many.pdf")
    g.savefig(path, bbox_inches="tight")

    # plt.figure()
    # sns.palplot(palette)
    # plt.savefig("palette.pdf")


@bacli.command
def folds(directory: str):
    # directory
    #   |- iteration 0
    #       |- train_fold_0.csv
    #       |- test_fold_0.csv
    #   |- iteration 1
    #       |- train_fold_1.csv
    #       |- test_fold_1.csv
    #   |- ...

    input_directory = Path(directory)

    test_fold_list = [
        i
        for i in input_directory.rglob("test_fold_*.csv")
        if not i.parent.name.startswith("_")
    ]

    df_list = []
    for i, fold in enumerate(test_fold_list):
        df = pd.read_csv(fold, sep=";")
        df["fold"] = i
        df_list.append(df)
    df = pd.concat(df_list).reset_index(drop=True)

    df = df.sort_values(by=["antigen.epitope", "y"]).reset_index(drop=True)

    # only plot image if there are fewer than 5 folds because
    # it takes too long otherwise and the image becomes too large
    if df["fold"].nunique() <= 5:
        plot_cv_folds(df, input_directory / "cv_folds.pdf")

    # save csv with test fold sample counts
    df.groupby(["fold", "antigen.epitope", "y"]).agg("size").groupby(
        level=0, group_keys=False
    ).apply(lambda x: x.sort_values(ascending=False).head(len(df))).to_csv(
        input_directory / "folds_pos_neg_split.csv"
    )

    df.groupby(["fold", "antigen.epitope"]).agg("size").groupby(
        level=0, group_keys=False
    ).apply(lambda x: x.sort_values(ascending=False).head(len(df))).to_csv(
        input_directory / "folds_pos_neg.csv"
    )

    df[df["y"] == 1].groupby(["fold", "antigen.epitope"]).agg("size").groupby(
        level=0, group_keys=False
    ).apply(
        lambda x: x.sort_values(ascending=False).head(len(df[df["y"] == 1]))
    ).to_csv(
        input_directory / "folds_pos_only.csv"
    )


@bacli.command
def metrics(directory: str, force: bool = False, y_lim_loss: float = None):
    # directory
    #   |- iteration 0
    #       |- metrics.csv
    #       |- roc.csv
    #       |- precision_recall.csv
    #       |- average_precision.csv
    #   |- iteration 1
    #       |- ...
    derive_metrics_all(directory, force=force)
    consolidate_all(directory, force=force)
    plot_all(directory, y_lim_loss=y_lim_loss)


@bacli.command
def compare(root_dir: str, force: bool = False, y_lim_loss: float = None):
    for directory in subdirs(root_dir):
        if os.path.basename(directory).startswith("_"):
            print(
                "Found directory:", directory, "which starts with underscore, skipping"
            )
            continue

        print(f"Consolidating metrics of {directory}")
        consolidate_all(directory, force=force)

    print(f"Concatenating metrics of root: {root_dir}")
    concatenate_all(root_dir, force=force)

    print("Plotting")
    plot_all(root_dir, y_lim_loss=y_lim_loss)


@bacli.command
def joinplot(directories: list):
    plot_combined([d for d in directories if not os.path.basename(d).startswith("_")])


@bacli.command
def evaluate_self_plots(
    root_dir: str,
    min_obs: int = 30,
    min_iterations: int = 3,
    grouped: bool = False,
    decoy: bool = False,
):

    # warn about grouped setting
    if grouped:
        print("Ignoring min_iterations since grouped=True.")
        min_iterations = 1

    # check input directory
    input_directory = root_dir
    # input_directory = Path(args.input)
    # assert input_directory.is_dir(), "Input is not a directory..."

    # find metrics_per_epitope.csv
    per_epitope_df = pd.read_csv(
        os.path.join(input_directory, "metrics_per_epitope.csv")
    )

    per_epitope_df["type"] = Path(input_directory).name

    # get number of testing data points
    per_epitope_df["n"] = per_epitope_df.pos_data + per_epitope_df.neg_data

    # minimum number of data points required to include
    per_epitope_df = per_epitope_df[per_epitope_df["n"] >= min_obs].reset_index(
        drop=True
    )

    roc_per_epitope(
        eval_df=per_epitope_df,
        output_path=os.path.join(
            input_directory,
            "roc_per_epitope-"
            + str(Path(root_dir).absolute().name).replace(" ", "-")
            + ".pdf",
        ),
        min_iterations=min_iterations,
        grouped=grouped,
        decoy=decoy,
    )

    plt.close("all")

    roc_train_corr(
        eval_df=per_epitope_df,
        output_path=os.path.join(
            input_directory,
            "roc_train_correlation-"
            + str(Path(root_dir).absolute().name).replace(" ", "-")
            + ".pdf",
        ),
    )

    plt.close("all")

    roc_min_dist_box(
        eval_df=per_epitope_df,
        output_path=os.path.join(
            input_directory,
            "roc_min_dist_boxplot-"
            + str(Path(root_dir).absolute().name).replace(" ", "-")
            + ".pdf",
        ),
    )

    roc_avg_dist_corr(
        eval_df=per_epitope_df,
        output_path=os.path.join(
            input_directory,
            "roc_mean_dist_correlation-"
            + str(Path(root_dir).absolute().name).replace(" ", "-")
            + ".pdf",
        ),
    )

    plt.close("all")

    roc_min_dist_corr(
        eval_df=per_epitope_df,
        output_path=os.path.join(
            input_directory,
            "roc_min_dist_correlation-"
            + str(Path(root_dir).absolute().name).replace(" ", "-")
            + ".pdf",
        ),
    )

    plt.close("all")

    roc_median_dist_corr(
        eval_df=per_epitope_df,
        output_path=os.path.join(
            input_directory,
            "roc_median_dist_correlation-"
            + str(Path(root_dir).absolute().name).replace(" ", "-")
            + ".pdf",
        ),
    )

    plt.close("all")


@bacli.command
def evaluate_self_comparison_plots(
    root_dir: str,
    min_obs: int = 30,
    min_iterations: int = 3,
    grouped: bool = False,
    decoy: bool = False,
):
    """Create per-epitope comparison plots for different models.

    Parameters
    ----------
    root_dir : str
        Path to directory containing a subdirectory for each model to be compared.
    min_obs : int, optional
        The minimum number of examples an epitope needs before being included in the comparison, by default 30
    min_iterations : int, optional
        The minimum number of iterations (per model) an epitope needs to appear in before being included in the comparison, by default 3
    grouped : bool, optional
        Whether the underlying data came from an epitope-grouped model, by default False
    decoy : bool, optional
        Whether there are decoy and normal models to be compared, in which case the min_iterations across models should be split between decoy and normal, by default False. Not required when comparing only decoy models. Skips wilcoxon test.
    """
    # warn about grouped setting and do not set a minimum iterations limit.
    if grouped:
        print("Ignoring min_iterations since grouped=True.")
        min_iterations = 1

    df_list = list()

    for directory in subdirs(root_dir):
        if os.path.basename(directory).startswith("_"):
            print(
                "Found directory:", directory, "which starts with underscore, skipping"
            )
            continue

        print(
            f"Retrieving metrics per epitope from evaluate_test_folds subdirectory of {directory}"
        )

        eval_list = [
            i
            for i in directory.rglob("*metrics_per_epitope*.csv")
            if not i.parent.name.startswith("_")
        ]

        assert (
            len(eval_list) == 1
        ), f"Found no or multiple evaluate_test_folds.csv files in {directory}, aborting..."

        df = pd.read_csv(eval_list[0])
        df["type"] = directory.name
        df_list.append(df)

    df_concat = pd.concat(df_list)

    # get number of testing data points
    df_concat["n"] = df_concat.pos_data + df_concat.neg_data

    # minimum number of data points required to include
    df_concat = df_concat[df_concat["n"] >= min_obs].reset_index(drop=True)

    # create output name
    output_path = os.path.join(
        root_dir,
        "roc_per_epitope-"
        + "min_obs"
        + str(min_obs)
        + "-min_iterations"
        + str(min_iterations)
        + "-"
        + str(Path(root_dir).absolute().name).replace(" ", "-")
        + ".pdf",
    )

    roc_per_epitope(
        eval_df=df_concat,
        output_path=output_path,
        min_iterations=min_iterations,
        comparison=True,
        grouped=grouped,
        decoy=decoy,
    )

    plt.close("all")
