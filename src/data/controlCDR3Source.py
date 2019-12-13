import pandas as pd

from src.config import PROJECT_ROOT
from src.data.data_source import DataSource

CONTROL_CDR3_PATH = PROJECT_ROOT / "data/raw/CDR3_control_sequences.tsv"


class ControlCDR3Source(DataSource):
    def __init__(self, filepath=CONTROL_CDR3_PATH):
        super().__init__()
        self.filepath = filepath
        print("Reading ref CDR3")
        self.data = pd.read_csv(self.filepath, sep="\t").drop_duplicates(
            subset="CDR3_beta"
        )
        print("Done reading")

    def __len__(self):
        return len(self.data)

    def __iter__(self, *args, **kwargs):

        # lengths = defaultdict(int)
        for index, row in self.data.iterrows():
            pep = row["CDR3_beta"]
            # lengths[len(pep)] += 1
            yield pep

        # print("lengths")
        # for k, v in sorted(lengths.items()):
        #     print(k, v)
