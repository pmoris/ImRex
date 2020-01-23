import pandas as pd

from src.config import PROJECT_ROOT
from src.data.data_source import DataSource


VDJDB_PATH = PROJECT_ROOT / "data/interim/vdjdb-human-trb.csv"


class VdjdbSource(DataSource):
    """Object holding VDJDB data.
    Contains a pandas DataFrame and dictionary with header names.
    Implements an __iter__ method, and consequently can
    be iterated through via a loop or list comprehension to yield
    the cdr3 and epitope sequences as a tuple, plus a 1.

    Inherits from DataSource object,
    which in turn inherits from Stream object.
    """

    def __init__(
        self,
        filepath=VDJDB_PATH,
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    ):
        super().__init__()
        self.filepath = filepath
        self.data = pd.read_csv(self.filepath, sep=";")
        self.headers = headers

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        for index, row in self.data.iterrows():
            pep1 = row[self.headers["cdr3_header"]]
            pep2 = row[self.headers["epitope_header"]]
            yield (pep1, pep2), 1
