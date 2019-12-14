""" Scenario for neural network. """
import src.bacli as bacli
from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parseFeatures, parseOperator
from src.config import PROJECT_ROOT
from src.data.controlCDR3Source import ControlCDR3Source
from src.data.vdjdbSource import VdjdbSource
from src.models.modelPadded import ModelPadded
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
    operator: str = "absdiff",  # can be: prod, diff, absdiff, layer or best
    early_stop=False,
    data_path=PROJECT_ROOT / "data/interim/vdjdb-human-no10x.csv",
    neg_ref: bool = False,
    stratified: bool = False,
    optimizer: str = "rmsprop",
    learning_rate: bool = False,
):

    # print function arguments
    print(locals())

    dataSource = VdjdbSource(filepath=data_path)

    featuresList = parseFeatures(features)
    operator = parseOperator(operator)
    featureBuilder = CombinedPeptideFeatureBuilder(featuresList, operator)

    if stratified and val_split is not None:
        raise RuntimeError("Can't test stratified without k folds")

    print("features:", featuresList)
    print("operator:", operator)
    print("neg_ref:", neg_ref)
    print("stratified:", stratified)

    inverseMap = InverseMap()

    pep1Range = (min1, max1)
    pep2Range = (min2, max2)

    trainer = Trainer(epochs, lookup=inverseMap, includeEarlyStop=early_stop)
    model = ModelPadded(
        max1,
        max2,
        nameSuffix=name,
        channels=featureBuilder.getNumberLayers(),
        optimizer=optimizer,
        include_lr=learning_rate,
    )

    if val_split is not None:
        train, val = Splitter(dataSource, ratio=val_split)
        if neg_ref:
            negativeSource = ControlCDR3Source()
            negTrain, negVal = Splitter(negativeSource, ratio=val_split)
            iterations = [((train, negTrain), (val, negVal))]
        else:
            iterations = [(train, val)]

    else:
        FoldSplitter = (
            EpitopeStratifiedFoldSplitter if stratified else RandomFoldSplitter
        )
        folds = FoldSplitter(dataSource, nrFolds)
        iterations = FoldIterator(folds)
        if neg_ref:
            negativeSource = ControlCDR3Source()
            neg_folds = FoldSplitter(negativeSource, nrFolds)
            neg_iterations = FoldIterator(neg_folds)
            iterations = [
                ((train, negTrain), (val, negVal))
                for (train, val), (negTrain, negVal) in zip(iterations, neg_iterations)
            ]

    for index, (train, val) in enumerate(iterations):
        negTrain, negVal = None, None
        print("Iteration:", index)
        print("batch size", batch_size)
        if neg_ref:
            train, negTrain = train
            val, negVal = val
            print("neg train set", len(negTrain))
            print("neg val set", len(negVal))

        print("train set", len(train))
        print("val set", len(val))

        trainStream = PaddedBatchGenerator(
            dataStream=train,
            featureBuilder=featureBuilder,
            negRatio=neg_ratio,
            batchSize=batch_size,
            pep1Range=pep1Range,
            pep2Range=pep2Range,
            negativeStream=negTrain,
        )
        valStream = PaddedBatchGenerator(
            dataStream=val,
            featureBuilder=featureBuilder,
            negRatio=neg_ratio,
            batchSize=batch_size,
            pep1Range=pep1Range,
            pep2Range=pep2Range,
            inverseMap=inverseMap,
            negativeStream=negVal,
        )

        trainer.train(model, trainStream, valStream, iteration=index)
