import os

from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import src.bacli as bacli
from src.bio.image import imageFromMatrix, imageFromMatrices
from src.bio.peptide_feature import (
    Charge,
    Hydrophilicity,
    Hydrophobicity,
    Mass,
    Surface,
    Flexibility,
    Transfer,
    TCRexBasicity,
    TCRexHelicity,
    AtchleyFactor1,
    AtchleyFactor2,
    AtchleyFactor3,
    AtchleyFactor4,
    AtchleyFactor5,
    parseFeatures,
    parseOperator,
    Polarity,
)
from src.bio.util import subdirs
from src.config import PROJECT_ROOT
from src.definitions.amino_acid_properties import AMINO_ACIDS
from src.metric import metric
from src.visualisation.plot import (
    consolidateAll,
    concatenateAll,
    plotAll,
    # palette,
    cmap,
    plotCombined,
)


# OUTPUT_DIR = "output/"
OUTPUT_DIR = PROJECT_ROOT / "reports/figures"
SCALE = 50

dependencies = {
    "mean_pred": metric.mean_pred,
    "AUC": metric.AUC,
    "balanced_accuracy": metric.balanced_accuracy,
}


@bacli.command
def test(modelFile: str):
    """ Debug function used to test commands out. """
    from keras.models import load_model

    model = load_model(modelFile, custom_objects=dependencies)
    print(model.layers[0].get_weights())


@bacli.command
def render(modelFile: str):
    """ Render dot representation of neural network. """
    from keras.models import load_model
    from keras.utils.vis_utils import plot_model

    model = load_model(modelFile, compile=False)
    plot_model(
        model,
        to_file=modelFile + ".pdf",
        show_shapes=True,
        show_layer_names=True,
        rankdir="TB",
    )


@bacli.command
def activations(
    modelFile: str,
    epitope: str = None,
    cdr3: str = None,
    layer: str = "conv2d_1",
    vsc: bool = False,
    features: str = "mass,hydrophob,charge,hydrophil",
    operator: str = "best",
):
    """ Display the activations of a model for a given input. """
    from keras.models import load_model
    from keract import get_activations, display_activations  # , display_heatmaps

    from src.processing.image_padding import ImagePadding
    from src.bio.feature_builder import CombinedPeptideFeatureBuilder

    if epitope is None or cdr3 is None:
        print("Supply epitope and cdr3")
        exit()

    model = load_model(modelFile, compile=False)
    if vsc:
        model.summary(line_length=80)
        model = model.get_layer("sequential_1")
    model.summary(line_length=80)
    model.compile(loss="mse", optimizer="adam")

    featuresList = parseFeatures(features)
    operator = parseOperator(operator)
    featureBuilder = CombinedPeptideFeatureBuilder(featuresList, operator)
    image = featureBuilder.generateFeature([cdr3, epitope])
    padded, _ = ImagePadding(None, 20, 13).transform((image, None), padValue=0)

    x = np.array((padded,))
    activations = get_activations(model, x, layer)
    display_activations(activations, cmap="gray", save=False)
    # display_heatmaps(activations, padded, save=False)


@bacli.command
def summary(modelFile: str):
    """ Print summary of neural network. """
    from keras.models import load_model

    model = load_model(modelFile, custom_objects=dependencies)
    model.summary(line_length=80)


@bacli.command
def peptide(
    epitope: str = None,
    cdr3: str = None,
    tensor: bool = False,
    cmyk: bool = False,
    operator: str = "best",
):
    """ Render peptide images. """
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for pep in [epitope, cdr3]:
        print(pep)

    mode = "CMYK" if cmyk else "RGB"

    matrices = list()
    layers = list()

    features = [Mass(), Hydrophobicity(), Charge()]
    if cmyk:
        features.append(Hydrophilicity())

    operator = parseOperator(operator)

    for index, feature in enumerate(features):
        matrix = feature.image_matrix(cdr3, epitope, operator=operator)
        print(matrix)
        matrices.append(matrix)
        if tensor:
            pixels = feature.image_tensor(cdr3, epitope)
        else:
            pixels = feature.image_matrix(cdr3, epitope, operator=operator)

        img = imageFromMatrix(pixels, mode=mode, index=index)
        layers.append((img, feature.name))

    img = imageFromMatrices(*matrices, mode=mode)
    layers.append((img, "Combined"))

    imgToPlot(layers, epitope, cdr3, "Image.pdf")


def imgToPlot(layers, epitope, cdr3, name):
    fig = plt.figure(figsize=(12, 4))
    axes = fig.subplots(1, 5, sharey=True)
    for ax in axes.flat:
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
    for ax in axes.flat:
        ax.label_outer()

    for i, (layer, title) in enumerate(layers):
        sub = axes[i]
        sub.set_title(title)
        rgb_im = layer.convert("RGB")
        pix = np.array(rgb_im)
        sub.imshow(pix, origin="lower")
        sub.grid(False)
        for spine in plt.gca().spines.values():
            spine.set_visible(False)

    plt.savefig(os.path.join(OUTPUT_DIR, name), bbox_inches="tight")


@bacli.command
def features():
    data = list()
    features = [
        Charge(),
        Hydrophobicity(),
        Hydrophilicity(),
        Polarity(),
        Mass(),
        Surface(),
        Flexibility(),
        Transfer(),
        TCRexBasicity(),
        TCRexHelicity(),
        AtchleyFactor1(),
        AtchleyFactor2(),
        AtchleyFactor3(),
        AtchleyFactor4(),
        AtchleyFactor5(),
    ]
    # features = [Charge(), Hydrophobicity(), Polarity(), Mass(), Hydrophilicity()]
    # features = [Charge(), Hydrophobicity(), Polarity(), Mass(), Hydrophilicity(), TCRexBasicity(), TCRexHelicity(), TCRexHydrophobicity(), TCRexMutationStability()]
    for f in features:
        print(f.name)
        values = list()
        for aa in AMINO_ACIDS:
            values.append(f.value(aa))
            print(f"{aa};{f.value(aa)}")
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
    path = os.path.join(OUTPUT_DIR, "features.pdf")
    g.savefig(path, bbox_inches="tight")

    # plt.figure()
    # sns.palplot(palette)
    # plt.savefig("palette.pdf")


@bacli.command
def metrics(directory: str, force: bool = False):
    # directory
    #   |- iteration 0
    #       |- metrics.csv
    #       |- roc.csv
    #       |- precision_recall.csv
    #       |- average_precision.csv
    #   |- iteration 1
    #       |- ...
    consolidateAll(directory, force=force)
    plotAll(directory)


@bacli.command
def compare(rootDirectory: str, force: bool = False):
    for directory in subdirs(rootDirectory):
        if os.path.basename(directory).startswith("_"):
            print(
                "Found directory:", directory, "which starts with underscore, skipping"
            )
            continue

        print(f"Consolidating metrics of {directory}")
        consolidateAll(directory, force=force)

    print(f"Concatenating metrics of root: {rootDirectory}")
    concatenateAll(rootDirectory, force=force)

    print("Plotting")
    plotAll(rootDirectory)


@bacli.command
def joinplot(directories: list):
    plotCombined([d for d in directories if not os.path.basename(d).startswith("_")])
