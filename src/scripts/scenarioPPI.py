""" Scenario for neural network. """
import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

import bacli

from data.ppiSource import PpiSource, SequencesMap
from models.modelPPIPadded import ModelPPIPadded
from bio.feature_builder import *
from processing.padded_batch_generator import PaddedBatchGenerator, PaddedBatchGenerator2
from processing.kfolds import FoldIterator, RandomFoldSplitter
from processing.splitter import Splitter
from neural.trainer import Trainer


bacli.setDescription(__doc__)


@bacli.command
def run(batch_size: int=16,
        val_split: float=None,
        epochs: int=40,
        neg_ratio: float=0.5,
        min: int=0,
        max: int=600,
        name: str="",
        nrFolds: int=3,
        features: str = "hydrophob,polarity,mass,hydrophil,charge",
        operator: str = "best",  # can be: prod, diff, layer or best
        early_stop=False,
        data_path='../data/PPI_positive.csv',
        negative_path='../data/PPI_negative.csv',
        sequences_path='../data/PPI_sequences.csv',
        swap=False,
        ):

    featuresList = parseFeatures(features)
    operator = parseOperator(operator)
    featureBuilder = CombinedPeptideFeatureBuilder(featuresList, operator)

    sequencesMap = SequencesMap(sequences_path)
    ppiSourcePos = PpiSource(data_path, sequencesMap, label=1)

    print("features:", featuresList)
    print("operator:", operator)
    print("swap:", swap)

    pepRange = (min, max)

    trainer = Trainer(epochs, includeEarlyStop=early_stop)
    model = ModelPPIPadded(max, max, nameSuffix=name, channels=featureBuilder.getNumberLayers())

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
            iterations = [((train, negTrain), (val, negVal)) for (train, val), (negTrain, negVal) in zip(iterations, neg_iterations)]

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

            trainStream = PaddedBatchGenerator2(train, negTrain, featureBuilder, neg_ratio, batch_size, pepRange, pepRange, swap=swap)
            valStream = PaddedBatchGenerator2(val, negVal, featureBuilder, neg_ratio, batch_size, pepRange, pepRange, swap=swap)

        else:
            print("train set", len(train))
            print("val set", len(val))

            trainStream = PaddedBatchGenerator(train, featureBuilder, neg_ratio, batch_size, pepRange, pepRange, cacheImages=False, swap=swap)
            valStream = PaddedBatchGenerator(val, featureBuilder, neg_ratio, batch_size, pepRange, pepRange, cacheImages=False, swap=swap)

        trainer.train(model, trainStream, valStream, iteration=index)


# if __name__ == "__main__":
#     run(val_split=0.2, batch_size=128, features="mass,charge,hydrophob", operator="best", data_path="../data/PPI_positive_SUBSET.csv")
#     exit()

