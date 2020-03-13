import os

from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import src.bacli as bacli
from src.bio.image import image_from_matrix, image_from_matrices
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
    parse_features,
    parse_operator,
    IsoelectricPoint,
)
from src.bio.util import subdirs
from src.config import PROJECT_ROOT
from src.definitions.amino_acid_properties import AMINO_ACIDS
from src.metric import metric
from src.visualisation.plot import (
    consolidate_all,
    concatenate_all,
    plot_all,
    # palette,
    cmap,
    plot_combined,
)


# OUTPUT_DIR = "output/"
OUTPUT_DIR = PROJECT_ROOT / "reports/figures"
SCALE = 50

dependencies = {
    "mean_pred": metric.mean_pred,
    "AUC": metric.auc,
    "balanced_accuracy": metric.balanced_accuracy,
}


@bacli.command
def test(model_file: str):
    """ Debug function used to test commands out. """
    from keras.models import load_model

    model = load_model(model_file, custom_objects=dependencies)
    print(model.layers[0].get_weights())


@bacli.command
def render(model_file: str):
    """ Render dot representation of neural network. """
    from keras.models import load_model
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
    from keras.models import load_model
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


def makeInput(epitope, cdr3):
    # TODO: implement
    #  generate input image for epitope and cdr3 sequence
    #  need to handle 'X' amino acid
    return None


@bacli.command
def cam(model_file: str, epitope, cdr3):
    from keras.models import load_model
    from vis.visualization import visualize_cam

    model = load_model(model_file, custom_objects=dependencies)

    # last layer is combination of previous one
    # originalLayer, title = getImage(epitope, cdr3, False, True, 'best')[-1]

    x = makeInput(epitope, cdr3)

    for modifier in [None, 'guided', 'relu']:
        cam = visualize_cam(model,
                            -1,
                            None,
                            x,
                            penultimate_layer_idx=None,
                            backprop_modifier=modifier,
                            grad_modifier=None
                            )
        camLayer = image_from_matrix(cam)

        if modifier is None:
            modifier = 'vanilla'

        layers = [
            # (originalLayer, "input", "CMYK"),
            (camLayer, modifier, "jet")
        ]

        img2plot(layers, epitope, cdr3, "CAM.pdf")


@bacli.command
def predict_variations(model_file: str, epitope, cdr3):
    from keras.models import load_model

    model = load_model(model_file, custom_objects=dependencies)

    samples = list()

    # first sample is original
    samples.append((epitope, cdr3))

    # insert permutations of epitope
    for i in range(len(epitope)):
        var_epitope = epitope.copy()
        var_epitope[i] = 'X'
        samples.append((var_epitope, cdr3))

    # insert permutations of cdr3
    for i in range(len(epitope)):
        var_cdr3 = cdr3.copy()
        var_cdr3[i] = 'X'
        samples.append((epitope, var_cdr3))

    # generate input images for each sample
    x = [makeInput(*sample) for sample in samples]

    # get scores of model
    predictions = model.predict_on_batch(x)

    # output results
    print("=============== All results ===============")
    for sample, prediction in zip(samples, predictions):
        print(sample[0], sample[1], prediction)

    print("=============== Statistics ===============")
    basePrediction = predictions[0]
    mostDifferent = sorted(zip(samples, predictions), key=lambda el: abs(el[1] - basePrediction), reverse=True)[5]

    print("Base prediction:", basePrediction)
    print("5 Most deviating:")
    for sample, prediction in mostDifferent:
        print("\t", sample[0], sample[1], prediction)


@bacli.command
def summary(model_file: str):
    """ Print summary of neural network. """
    from keras.models import load_model

    model = load_model(model_file, custom_objects=dependencies)
    model.summary(line_length=80)


def getImage(epitope: str = None,
             cdr3: str = None,
             tensor: bool = False,
             cmyk: bool = False,
             operator: str = "best",):

    for pep in [epitope, cdr3]:
        print(pep)

    mode = "CMYK" if cmyk else "RGB"

    matrices = list()
    layers = list()

    features = [Mass(), Hydrophobicity(), Charge()]
    if cmyk:
        features.append(Hydrophilicity())

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
        layers.append((img, feature.name, mode))

    img = image_from_matrices(*matrices, mode=mode)
    layers.append((img, "Combined", mode))
    return layers


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

    layers = getImage(epitope, cdr3, tensor, cmyk, operator)

    img2plot(layers, epitope, cdr3, "amino-acid-map.pdf")


def img2plot(layers, epitope, cdr3, name):
    fig = plt.figure(figsize=(12, 4))
    axes = fig.subplots(1, len(layers), sharey=True)
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

    for i, (layer, title, cmap) in enumerate(layers):
        sub = axes[i]
        sub.set_title(title)
        rgb_im = layer.convert("RGB")
        pix = np.array(rgb_im)
        sub.imshow(pix, origin="lower", cmap=cmap)
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
def metrics(directory: str, force: bool = False):
    # directory
    #   |- iteration 0
    #       |- metrics.csv
    #       |- roc.csv
    #       |- precision_recall.csv
    #       |- average_precision.csv
    #   |- iteration 1
    #       |- ...
    consolidate_all(directory, force=force)
    plot_all(directory)


@bacli.command
def compare(root_dir: str, force: bool = False):
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
    plot_all(root_dir)


@bacli.command
def joinplot(directories: list):
    plot_combined([d for d in directories if not os.path.basename(d).startswith("_")])
