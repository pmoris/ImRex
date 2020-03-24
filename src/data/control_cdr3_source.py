import logging
from typing import Optional

import pandas as pd

from src.config import PROJECT_ROOT
from src.data.data_source import DataSource

CONTROL_CDR3_PATH = PROJECT_ROOT / "data/raw/CDR3_control_sequences.tsv"
CONTROL_CDR3_SEQ_COLUMN = "CDR3_beta"


def restrict_length(df: pd.DataFrame, min_length: int, max_length: int):
    df = df.loc[
        (df[CONTROL_CDR3_SEQ_COLUMN].str.len() >= min_length)
        & (df[CONTROL_CDR3_SEQ_COLUMN].str.len() <= max_length)
    ]
    return df


class ControlCDR3Source(DataSource):
    def __init__(
        self,
        filepath=CONTROL_CDR3_PATH,
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
        headers={"cdr3": CONTROL_CDR3_SEQ_COLUMN},
    ):
        super().__init__()
        self.filepath = filepath
        self.headers = headers

        logger = logging.getLogger(__name__)

        logger.info(f"Read reference CDR3 sequences from '{filepath}'")

        # set general min and max values if not provided
        if not min_length:
            logger.info(f"Minimum sequence length was not provided, defaulting to 0.")
            min_length = 0
        if not max_length:
            logger.info(f"Maximum sequence length was not provided, defaulting to 100.")
            max_length = 100

        # read data
        self.data = (
            pd.read_csv(self.filepath, sep="\t")
            .pipe(restrict_length, min_length=min_length, max_length=max_length)
            .drop_duplicates(subset=CONTROL_CDR3_SEQ_COLUMN)
        )

        logger.info(
            f"Retrieved {self.data.shape[0]} unique CDR3 sequences with a length between {min_length} and {max_length}"
        )

    def __len__(self):
        return len(self.data)

    def __iter__(self, *args, **kwargs):
        for index, row in self.data.iterrows():
            pep = row[CONTROL_CDR3_SEQ_COLUMN]
            yield pep
