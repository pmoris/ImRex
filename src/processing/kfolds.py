import random
from collections import defaultdict

from src.processing.data_stream import DataStream


def epitope_stratified_fold_splitter(data_source, n_folds):
    epitope_sets = defaultdict(list)
    for (cdr3, epitope), label in data_source:
        epitope_sets[epitope].append((cdr3, label))

    folds = [list() for _ in range(n_folds)]
    epitope_counts = [0] * n_folds
    for epitope, values in sorted(
        epitope_sets.items(), key=lambda x: len(x[1]), reverse=True
    ):

        # add epitope to fold with least entries
        min_amount = len(min(folds, key=lambda x: len(x)))
        candidates = [
            (index, f, epitope_count)
            for index, (f, epitope_count) in enumerate(zip(folds, epitope_counts))
            if len(f) == min_amount
        ]

        # select fold with least distinct epitopes if tie
        index, fold, epitope_count = min(candidates, key=lambda x: x[2])
        assert epitope_count == epitope_counts[index]

        fold += [((cdr3, epitope), label) for cdr3, label in values]
        epitope_counts[index] += 1

    print("Folds:")
    for f, epitopes in zip(folds, epitope_counts):
        print("size:", len(f), "epitopes:", epitopes)

    for epitopes in epitope_counts:
        if epitopes < 2:
            raise RuntimeError(
                "Can't validate fold with one epitope. Need at least two for the generation of negative samples."
            )

    return folds


def random_fold_splitter(data_source, n_folds, shuffle=True):
    data = list(data_source)
    if shuffle:
        random.shuffle(data)
    folds = [data[i::n_folds] for i in range(n_folds)]

    print("Folds:")
    for f in folds:
        print("size:", len(f))

    return folds


def fold_iterator(folds):
    for leave_out in range(len(folds)):
        train = flatten(folds[:leave_out] + folds[leave_out + 1 :])
        val = folds[leave_out]
        yield DataStream(train), DataStream(val)


def flatten(l):
    return [item for sublist in l for item in sublist]
