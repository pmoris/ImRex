from functools import lru_cache

import pandas as pd

from src.data.data_source import DataSource


POS_PATH = "../data/raw/ppi/PPI_positive.csv"
# NEG_PATH = "../data/raw/ppi/PPI_negative.csv"
DICT_PATH = "../data/raw/ppi/PPI_sequences.csv"


@lru_cache()
class SequencesMap(object):
    def __init__(self, dict_path=DICT_PATH):
        data = pd.read_csv(dict_path, sep=";")
        self.map = {row[0]: row[1] for index, row in data.iterrows()}

    def get(self, key):
        value = self.map[key]
        return value


class PpiSource(DataSource):
    ID_A = "Unique identifier for interactor A"
    ID_B = "Unique identifier for interactor B"

    def __init__(self, filepath, sequences_map=SequencesMap(), label=None):
        super().__init__()
        self.filepath = filepath
        self.sequences_map = sequences_map
        self.data = pd.read_csv(self.filepath, sep=";")
        self.label = label

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        for index, row in self.data.iterrows():
            id_a = row[PpiSource.ID_A]
            id_b = row[PpiSource.ID_B]
            pep1 = self.sequences_map.get(id_a)
            pep2 = self.sequences_map.get(id_b)

            if self.label is not None:
                yield (pep1, pep2), self.label
            else:
                yield (pep1, pep2)


pos_ppi_source = PpiSource(POS_PATH, SequencesMap(), label=1)
# neg_ppi_source = PpiSource(NEG_PATH, SequencesMap())
