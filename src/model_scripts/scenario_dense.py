""" Scenario for neural network. """
import src.bacli as bacli

from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.models.model_dense import ModelDense
from src.neural.trainer import Trainer
from src.processing.inverse_map import InverseMap
from src.processing.kfolds import (
    epitope_stratified_fold_splitter,
    fold_iterator,
    random_fold_splitter,
)
from src.processing.padded_batch_generator import padded_batch_generator
from src.processing.splitter import splitter

bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 128,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min1: int = 10,
    max1: int = 20,
    min2: int = 8,
    max2: int = 13,
    name: str = "",
    n_folds: int = 5,
    features: str = "hydrophob,isoelectric,mass,hydrophil,charge",
    operator: str = "best",  # can be: prod,diff,layer or best
    early_stop=False,
    data_path=PROJECT_ROOT / "data/interim/vdjdb-human-no10x.csv",
    stratified: bool = False,
):

    data_source = VdjdbSource(filepath=data_path)

    features_list = parse_features(features)
    operator = parse_operator(operator)
    feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)

    print("features:", features_list)
    print("operator:", operator)

    inverse_map = InverseMap()

    pep1_range = (min1, max1)
    pep2_range = (min2, max2)

    trainer = Trainer(epochs, lookup=inverse_map, include_early_stop=early_stop)
    model = ModelDense(
        max1, max2, nameSuffix=name, channels=feature_builder.get_number_layers()
    )

    if val_split is not None:
        train, val = splitter(data_source, ratio=val_split)
        iterations = [(train, val)]
    else:
        fold_splitter = (
            epitope_stratified_fold_splitter if stratified else random_fold_splitter
        )
        folds = fold_splitter(data_source, n_folds)
        iterations = fold_iterator(folds)
    for index, (train, val) in enumerate(iterations):
        print("Iteration:", index)
        print("train set", len(train))
        print("val set", len(val))
        print("batch size", batch_size)

        train_stream = padded_batch_generator(
            train, feature_builder, neg_ratio, batch_size, pep1_range, pep2_range
        )
        val_stream = padded_batch_generator(
            val,
            feature_builder,
            neg_ratio,
            batch_size,
            pep1_range,
            pep2_range,
            inverse_map=inverse_map,
        )

        trainer.train(model, train_stream, val_stream, iteration=index)
