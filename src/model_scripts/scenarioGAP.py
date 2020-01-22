""" Scenario for neural network. """
import src.bacli as bacli
from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parseFeatures, parseOperator
from src.config import PROJECT_ROOT
from src.data.vdjdbSource import VdjdbSource
from src.models.modelGAP import ModelGAP
from src.neural.trainer import Trainer
from src.processing.grouped_batch_generator import (
    GroupedBatchGenerator,
    # GroupedBatchGenerator2,
)
from src.processing.kfolds import (
    EpitopeStratifiedFoldSplitter,
    FoldIterator,
    RandomFoldSplitter,
)
from src.processing.splitter import Splitter

# from src.data.controlCDR3Source import ControlCDR3Source


bacli.setDescription(__doc__)


@bacli.command
def run(
    batch_size: int = 128,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min_group: int = 32,
    name: str = "",
    nrFolds: int = 5,
    features: str = "hydrophob,isoelectric,mass,hydrophil,charge",
    operator: str = "best",  # can be: prod,diff,layer or best
    early_stop=False,
    data_path=PROJECT_ROOT / "data/interim/vdjdb-human-no10x.csv",
    stratified: bool = False,
):

    dataSource = VdjdbSource(filepath=data_path)

    # negativeSource = ControlCDR3Source()
    # negTrain, negVal = Splitter(negativeSource, ratio=val_split)

    featuresList = parseFeatures(features)
    operator = parseOperator(operator)
    featureBuilder = CombinedPeptideFeatureBuilder(featuresList, operator)

    print("features:", featuresList)
    print("operator:", operator)

    trainer = Trainer(epochs, includeEarlyStop=early_stop)
    model = ModelGAP(nameSuffix=name, channels=featureBuilder.getNumberLayers())

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

        trainStream = GroupedBatchGenerator(
            train, featureBuilder, neg_ratio, batch_size, min_group
        )
        valStream = GroupedBatchGenerator(
            val, featureBuilder, neg_ratio, batch_size, min_group
        )

        # trainStream = PaddedBatchGenerator(train, featureBuilder, neg_ratio, batch_size, pep1Range, pep2Range)
        # valStream = PaddedBatchGenerator(val, featureBuilder, neg_ratio, batch_size, pep1Range, pep2Range, inverseMap=inverseMap)

        trainer.train(model, trainStream, valStream, iteration=index)
