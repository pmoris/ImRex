from collections import defaultdict
import random
from .data_stream import DataStream


def EpitopeStratifiedFoldSplitter(dataSource, nrFolds):
    epitopeSets = defaultdict(list)
    for (cdr3, epitope), label in dataSource:
        epitopeSets[epitope].append((cdr3, label))

    folds = [list() for _ in range(nrFolds)]
    epitopeCounts = [0] * nrFolds
    for epitope, values in sorted(
        epitopeSets.items(), key=lambda x: len(x[1]), reverse=True
    ):

        # add epitope to fold with least entries
        minAmount = len(min(folds, key=lambda x: len(x)))
        candidates = [
            (index, f, epitopeCount)
            for index, (f, epitopeCount) in enumerate(zip(folds, epitopeCounts))
            if len(f) == minAmount
        ]

        # select fold with least distinct epitopes if tie
        index, fold, epitopeCount = min(candidates, key=lambda x: x[2])
        assert epitopeCount == epitopeCounts[index]

        fold += [((cdr3, epitope), label) for cdr3, label in values]
        epitopeCounts[index] += 1

    print("Folds:")
    for f, epitopes in zip(folds, epitopeCounts):
        print("size:", len(f), "epitopes:", epitopes)

    for epitopes in epitopeCounts:
        if epitopes < 2:
            raise RuntimeError(
                "Can't validate fold with one epitope. Need at least two for the generation of negative samples."
            )

    return folds


def RandomFoldSplitter(dataSource, nrFolds, shuffle=True):
    data = list(dataSource)
    if shuffle:
        random.shuffle(data)
    folds = [data[i::nrFolds] for i in range(nrFolds)]

    print("Folds:")
    for f in folds:
        print("size:", len(f))

    return folds


def FoldIterator(folds):
    for leaveOut in range(len(folds)):
        train = flatten(folds[:leaveOut] + folds[leaveOut + 1 :])
        val = folds[leaveOut]
        yield DataStream(train), DataStream(val)


def flatten(l):
    return [item for sublist in l for item in sublist]
