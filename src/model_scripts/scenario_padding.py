""" Scenario for neural network. """
import src.bacli as bacli
from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.control_cdr3_source import ControlCDR3Source
from src.data.vdjdb_source import VdjdbSource
from src.models.model_padded import ModelPadded
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
    features: str = "hydrophob,isoelectric,mass,hydrophil,charge",  # can be any str listed in peptide_feature.featuresMap
    operator: str = "absdiff",  # can be: prod, diff, absdiff, layer or best
    early_stop=False,
    data_path=PROJECT_ROOT / "data/interim/vdjdb-human-no10x.csv",
    neg_ref: bool = False,
    stratified: bool = False,
    optimizer: str = "rmsprop",
    learning_rate: bool = False,
    dense_activation: str = "tanh",
):

    # print function arguments that were used for logging purposes
    print(locals())

    data_source = VdjdbSource(
        filepath=data_path,
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )

    features_list = parse_features(features)
    operator = parse_operator(operator)
    feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)

    if stratified and val_split is not None:
        raise RuntimeError("Can't test stratified without k folds")

    print("features:", features_list)
    print("operator:", operator)
    print("neg_ref:", neg_ref)
    print("stratified:", stratified)

    inverse_map = InverseMap()

    pep1_range = (min1, max1)
    pep2_range = (min2, max2)

    trainer = Trainer(epochs, lookup=inverse_map, include_early_stop=early_stop)
    model = ModelPadded(
        max1,
        max2,
        nameSuffix=name,
        channels=feature_builder.get_number_layers(),
        optimizer=optimizer,
        include_lr=learning_rate,
    )

    if val_split is not None:
        train, val = splitter(data_source, ratio=val_split)
        if neg_ref:
            negative_source = ControlCDR3Source()
            neg_train, neg_val = splitter(negative_source, ratio=val_split)
            iterations = [((train, neg_train), (val, neg_val))]
        else:
            iterations = [(train, val)]

    else:
        fold_splitter = (
            epitope_stratified_fold_splitter if stratified else random_fold_splitter
        )
        folds = fold_splitter(data_source, n_folds)
        iterations = fold_iterator(folds)
        if neg_ref:
            negative_source = ControlCDR3Source()
            neg_folds = fold_splitter(negative_source, n_folds)
            neg_iterations = fold_iterator(neg_folds)
            iterations = [
                ((train, neg_train), (val, neg_val))
                for (train, val), (neg_train, neg_val) in zip(
                    iterations, neg_iterations
                )
            ]

    for index, (train, val) in enumerate(iterations):
        neg_train, neg_val = None, None
        print("Iteration:", index)
        print("batch size", batch_size)
        if neg_ref:
            train, neg_train = train
            val, neg_val = val
            print("neg train set", len(neg_train))
            print("neg val set", len(neg_val))

        print("train set", len(train))
        print("val set", len(val))

        train_stream = padded_batch_generator(
            data_stream=train,
            feature_builder=feature_builder,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            pep1_range=pep1_range,
            pep2_range=pep2_range,
            negative_stream=neg_train,
        )
        val_stream = padded_batch_generator(
            data_stream=val,
            feature_builder=feature_builder,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            pep1_range=pep1_range,
            pep2_range=pep2_range,
            inverse_map=inverse_map,
            negative_stream=neg_val,
        )

        trainer.train(model, train_stream, val_stream, iteration=index)
