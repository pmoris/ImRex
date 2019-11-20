# from collections import defaultdict
from functools import lru_cache

import pandas as pd

from src.data.data_source import DataSource


POS_PATH = "../data/PPI_positive.csv"
# NEG_PATH = "../data/PPI_negative.csv"

DICT_PATH = "../data/PPI_sequences.csv"


@lru_cache()
class SequencesMap(object):
    def __init__(self, dictPath=DICT_PATH):
        data = pd.read_csv(dictPath, sep=";")
        self.map = {row[0]: row[1] for index, row in data.iterrows()}

    def get(self, key):
        value = self.map[key]
        return value


class PpiSource(DataSource):
    ID_A = "Unique identifier for interactor A"
    ID_B = "Unique identifier for interactor B"

    def __init__(self, filepath, sequencesMap=SequencesMap(), label=None):
        super().__init__()
        self.filepath = filepath
        self.sequencesMap = sequencesMap
        self.data = pd.read_csv(self.filepath, sep=";")
        self.label = label

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        # lengths = defaultdict(int)

        for index, row in self.data.iterrows():
            idA = row[PpiSource.ID_A]
            idB = row[PpiSource.ID_B]
            pep1 = self.sequencesMap.get(idA)
            pep2 = self.sequencesMap.get(idB)

            # lengths[len(pep1)] += 1
            # lengths[len(pep2)] += 1
            if self.label is not None:
                yield (pep1, pep2), self.label
            else:
                yield (pep1, pep2)

        # print("lengths")
        # for k, v in sorted(lengths.items()):
        #     print(k, v)


posPpiSource = PpiSource(POS_PATH, SequencesMap(), label=1)
# negPpiSource = PpiSource(NEG_PATH, SequencesMap())
