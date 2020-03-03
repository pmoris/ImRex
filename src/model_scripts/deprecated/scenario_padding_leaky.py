""" Scenario for neural network. """
import src.bacli as bacli
from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.control_cdr3_source import ControlCDR3Source
from src.data.vdjdb_source import VdjdbSource
from src.models.model_padded_leaky import ModelPaddedLeaky
from src.neural.trainer import Trainer
from src.processing.inverse_map import InverseMap
from src.processing.cv_folds import cv_splitter
from src.processing.padded_batch_generator import padded_batch_generator
from src.processing.splitter import splitter

bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 128,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min_length_cdr3: int = 10,
    max_length_cdr3: int = 20,
    min_length_epitope: int = 8,
    max_length_epitope: int = 13,
    name: str = "",
    n_folds: int = 5,
    features: str = "hydrophob,isoelectric,mass,hydrophil,charge",
    operator: str = "absdiff",  # can be: prod, diff, absdiff, layer or best
    early_stop=False,
    data_path=PROJECT_ROOT
    / "data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x.csv",
    neg_ref: bool = False,
    stratified: bool = False,
):

    # print function arguments
    print(locals())

    data_source = VdjdbSource(filepath=data_path)

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

    cdr3_range = (min_length_cdr3, max_length_cdr3)
    epitope_range = (min_length_epitope, max_length_epitope)

    trainer = Trainer(epochs, lookup=inverse_map, include_early_stop=early_stop)
    model = ModelPaddedLeaky(
        max_length_cdr3,
        max_length_epitope,
        name_suffix=name,
        channels=feature_builder.get_number_layers(),
    )

    if val_split is not None:
        train, val = splitter(data_source, test_size=val_split)
        if neg_ref:
            negative_source = ControlCDR3Source(
                min_length=min_length_cdr3, max_length=max_length_cdr3
            )
            neg_train, neg_val = splitter(negative_source, test_size=val_split)
            iterations = [((train, neg_train), (val, neg_val))]
        else:
            iterations = [(train, val)]

    else:
        fold_splitter = epitope_grouped_fold_splitter if stratified else k_fold_splitter
        folds = fold_splitter(data_source, n_folds)
        iterations = fold_iterator(folds)
        if neg_ref:
            negative_source = ControlCDR3Source(
                min_length=min_length_cdr3, max_length=max_length_cdr3
            )
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
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
            negative_stream=neg_train,
        )
        val_stream = padded_batch_generator(
            data_stream=val,
            feature_builder=feature_builder,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
            inverse_map=inverse_map,
            negative_stream=neg_val,
        )

        trainer.train(model, train_stream, val_stream, iteration=index)
