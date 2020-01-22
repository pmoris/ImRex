""" Scenario for neural network. """
import src.bacli as bacli
from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.models.model_gap import ModelGAP
from src.neural.trainer import Trainer
from src.processing.grouped_batch_generator import (
    grouped_batch_generator,
    # GroupedBatchGenerator2,
)
from src.processing.kfolds import (
    epitope_stratified_fold_splitter,
    fold_iterator,
    random_fold_splitter,
)
from src.processing.splitter import splitter

# from src.data.controlCDR3Source import ControlCDR3Source


bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 128,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min_group: int = 32,
    name: str = "",
    n_folds: int = 5,
    features: str = "hydrophob,isoelectric,mass,hydrophil,charge",
    operator: str = "best",  # can be: prod,diff,layer or best
    early_stop=False,
    data_path=PROJECT_ROOT / "data/interim/vdjdb-human-no10x.csv",
    stratified: bool = False,
):

    data_source = VdjdbSource(filepath=data_path)

    # negativeSource = ControlCDR3Source()
    # negTrain, negVal = Splitter(negativeSource, ratio=val_split)

    features_list = parse_features(features)
    operator = parse_operator(operator)
    feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)

    print("features:", features_list)
    print("operator:", operator)

    trainer = Trainer(epochs, include_early_stop=early_stop)
    model = ModelGAP(nameSuffix=name, channels=feature_builder.get_number_layers())

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

        train_stream = grouped_batch_generator(
            data_stream=train,
            feature_builder=feature_builder,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            min_amount=min_group,
        )
        val_stream = grouped_batch_generator(
            data_stream=val,
            feature_builder=feature_builder,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            min_amount=min_group,
        )

        # train_stream = PaddedBatchGenerator(train, feature_builder, neg_ratio, batch_size, pep1Range, pep2Range)
        # val_stream = PaddedBatchGenerator(val, feature_builder, neg_ratio, batch_size, pep1Range, pep2Range, inverseMap=inverseMap)

        trainer.train(model, train_stream, val_stream, iteration=index)
