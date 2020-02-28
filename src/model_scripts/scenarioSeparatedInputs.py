""" Scenario for neural network. """
import src.bacli as bacli
from src.config import PROJECT_ROOT
from src.data.controlCDR3Source import ControlCDR3Source
from src.data.vdjdbSource import VdjdbSource
from src.models.modelSeparatedInputs import ModelSeparatedInputs
from src.neural.trainer import Trainer
from src.processing.kfolds import (
    EpitopeStratifiedFoldSplitter,
    FoldIterator,
    RandomFoldSplitter,
)
from src.processing.net_tcr_batch_generator import NetTCRBatchGenerator
from src.processing.splitter import Splitter

bacli.setDescription(__doc__)


@bacli.command
def run(
    batch_size: int = 128,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min_group: int = 32,
    # min1: int = 10,
    # max1: int = 20,
    # min2: int = 8,
    # max2: int = 13,
    name: str = "",
    nrFolds: int = 5,
    early_stop=False,
    data_path=PROJECT_ROOT / "data/interim/vdjdb-human-trb.csv",
    # neg_ref: bool = False,
    stratified: bool = False,
    optimizer: str = "rmsprop",
    learning_rate: bool = False,
):

    dataSource = VdjdbSource(filepath=data_path)

    if stratified and val_split is not None:
        raise RuntimeError("Can't test stratified without k folds")

    # print("neg_ref:", neg_ref)
    print("stratified:", stratified)

    # pep1Range = (min1, max1)
    # pep2Range = (min2, max2)

    trainer = Trainer(epochs, includeEarlyStop=early_stop)
    model = ModelSeparatedInputs(
        optimizer=optimizer, include_lr=learning_rate, nameSuffix=name
    )

    if val_split is not None:
        train, val = Splitter(dataSource, ratio=val_split)
        iterations = [(train, val)]
        # if neg_ref:
        #     negativeSource = ControlCDR3Source()
        #     negTrain, negVal = Splitter(negativeSource, ratio=val_split)
        #     iterations = [((train, negTrain), (val, negVal))]
        # else:
        #     iterations = [(train, val)]
    else:
        FoldSplitter = (
            EpitopeStratifiedFoldSplitter if stratified else RandomFoldSplitter
        )
        folds = FoldSplitter(dataSource, nrFolds)
        iterations = FoldIterator(folds)
        # if neg_ref:
        #     negativeSource = ControlCDR3Source()
        #     neg_folds = FoldSplitter(negativeSource, nrFolds)
        #     neg_iterations = FoldIterator(neg_folds)
        #     iterations = [
        #         ((train, negTrain), (val, negVal))
        #         for (train, val), (negTrain, negVal) in zip(iterations, neg_iterations)
        #     ]

    for index, (train, val) in enumerate(iterations):
        print("Iteration:", index)
        print("train set", len(train))
        print("val set", len(val))
        print("batch size", batch_size)

        # negTrain, negVal = None, None
        # print("Iteration:", index)
        # print("batch size", batch_size)
        # if neg_ref:
        #     train, negTrain = train
        #     val, negVal = val
        #     print("neg train set", len(negTrain))
        #     print("neg val set", len(negVal))

        # print("train set", len(train))
        # print("val set", len(val))

        trainStream = NetTCRBatchGenerator(
            dataStream=train,
            negRatio=neg_ratio,
            batchSize=batch_size,
            # pep1Range=pep1Range,
            # pep2Range=pep2Range,
            minAmount=min_group,
            # negativeStream=negTrain,
        )
        valStream = NetTCRBatchGenerator(
            dataStream=val,
            negRatio=neg_ratio,
            batchSize=batch_size,
            # pep1Range=pep1Range,
            # pep2Range=pep2Range,
            minAmount=min_group,
            # negativeStream=negVal,
        )

        trainer.train(model, trainStream, valStream, iteration=index)
