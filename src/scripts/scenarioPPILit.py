""" Scenario for neural network. """
import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

import bacli

from data.ppiSource import PpiSource, SequencesMap
from models.modelPPILit import ModelPPILit
from processing.ppi_lit_generator import PPILitGenerator, PPILitGenerator2
from processing.kfolds import FoldIterator, RandomFoldSplitter
from processing.splitter import Splitter
from neural.trainer import Trainer


bacli.setDescription(__doc__)


@bacli.command
def run(
    batch_size: int = 512,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min: int = 0,
    max: int = 1200,
    name: str = "",
    nrFolds: int = 3,
    early_stop=False,
    data_path="../data/PPI_positive.csv",
    negative_path="../data/PPI_negative.csv",
    sequences_path="../data/PPI_sequences.csv",
    swap=False,
):

    sequencesMap = SequencesMap(sequences_path)
    ppiSourcePos = PpiSource(data_path, sequencesMap, label=1)
    pepRange = (min, max)

    print("swap:", swap)

    trainer = Trainer(epochs, includeEarlyStop=early_stop)
    model = ModelPPILit(max, max, nameSuffix=name)

    if val_split is not None:
        train, val = Splitter(ppiSourcePos, ratio=val_split)
        if negative_path:
            ppiSourceNeg = PpiSource(negative_path, sequencesMap, label=0)
            negTrain, negVal = Splitter(ppiSourceNeg, ratio=val_split)
            iterations = [((train, negTrain), (val, negVal))]
        else:
            iterations = [(train, val)]
    else:
        folds = RandomFoldSplitter(ppiSourcePos, nrFolds)
        iterations = FoldIterator(folds)
        if negative_path:
            ppiSourceNeg = PpiSource(negative_path, sequencesMap, label=0)
            neg_folds = RandomFoldSplitter(ppiSourceNeg, nrFolds)
            neg_iterations = FoldIterator(neg_folds)
            iterations = [
                ((train, negTrain), (val, negVal))
                for (train, val), (negTrain, negVal) in zip(iterations, neg_iterations)
            ]

    for index, (train, val) in enumerate(iterations):
        print("Iteration:", index)
        print("batch size", batch_size)
        if negative_path:
            train, negTrain = train
            val, negVal = val
            print("train set", len(train))
            print("val set", len(val))

            print("neg train set", len(negTrain))
            print("neg val set", len(negVal))

            trainStream = PPILitGenerator2(
                train, negTrain, neg_ratio, batch_size, pepRange, pepRange, swap=swap
            )
            valStream = PPILitGenerator2(
                val, negVal, neg_ratio, batch_size, pepRange, pepRange, swap=swap
            )

        else:
            print("train set", len(train))
            print("val set", len(val))

            trainStream = PPILitGenerator(
                train, neg_ratio, batch_size, pepRange, pepRange
            )
            valStream = PPILitGenerator(val, neg_ratio, batch_size, pepRange, pepRange)

        trainer.train(model, trainStream, valStream, iteration=index)


# if __name__ == "__main__":
#     run(val_split=0.2, batch_size=4, data_path="../data/PPI_positive_SUBSET.csv")
#     exit()
