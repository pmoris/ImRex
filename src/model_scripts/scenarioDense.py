""" Scenario for neural network. """
import src.bacli as bacli

from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parseFeatures, parseOperator
from src.data.vdjdbSource import VdjdbSource
from src.models.modelDense import ModelDense
from src.neural.trainer import Trainer
from src.processing.inverse_map import InverseMap
from src.processing.kfolds import (
    EpitopeStratifiedFoldSplitter,
    FoldIterator,
    RandomFoldSplitter,
)
from src.processing.padded_batch_generator import PaddedBatchGenerator
from src.processing.splitter import Splitter

bacli.setDescription(__doc__)


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
    nrFolds: int = 5,
    features: str = "hydrophob,polarity,mass,hydrophil,charge",
    operator: str = "best",  # can be: prod,diff,layer or best
    early_stop=False,
    data_path="../data/vdjdb_TRB.csv",
    stratified: bool = False,
):

    dataSource = VdjdbSource(filepath=data_path)

    featuresList = parseFeatures(features)
    operator = parseOperator(operator)
    featureBuilder = CombinedPeptideFeatureBuilder(featuresList, operator)

    print("features:", featuresList)
    print("operator:", operator)

    inverseMap = InverseMap()

    pep1Range = (min1, max1)
    pep2Range = (min2, max2)

    trainer = Trainer(epochs, lookup=inverseMap, includeEarlyStop=early_stop)
    model = ModelDense(
        max1, max2, nameSuffix=name, channels=featureBuilder.getNumberLayers()
    )

    if val_split is not None:
        train, val = Splitter(dataSource, ratio=val_split)
        iterations = [(train, val)]
    else:
        FoldSplitter = (
            EpitopeStratifiedFoldSplitter if stratified else RandomFoldSplitter
        )
        folds = FoldSplitter(dataSource, nrFolds)
        iterations = FoldIterator(folds)
    for index, (train, val) in enumerate(iterations):
        print("Iteration:", index)
        print("train set", len(train))
        print("val set", len(val))
        print("batch size", batch_size)

        trainStream = PaddedBatchGenerator(
            train, featureBuilder, neg_ratio, batch_size, pep1Range, pep2Range
        )
        valStream = PaddedBatchGenerator(
            val,
            featureBuilder,
            neg_ratio,
            batch_size,
            pep1Range,
            pep2Range,
            inverseMap=inverseMap,
        )

        trainer.train(model, trainStream, valStream, iteration=index)