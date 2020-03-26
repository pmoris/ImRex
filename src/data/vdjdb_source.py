import pandas as pd

from src.config import PROJECT_ROOT
from src.data.data_source import DataSource
from src.processing.negative_sampler import sample_pairs


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
        self.data["y"] = 1

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        for index, row in self.data.iterrows():
            pep1 = row[self.headers["cdr3_header"]]
            pep2 = row[self.headers["epitope_header"]]
            label = row["y"]
            yield (pep1, pep2), label

    def generate_negatives_from_ref(self, negative_source: DataSource):
        """ Generate negative CDR3 epitope sequence pairs by drawing from a negative CDR3 reference set.

        Every epitope in the positive set is matched with a random CDR3 in order to keep the epitope distribution equal between the two classes.
        """
        # sample required number of CDR3 sequences
        amount = self.data.shape[0]
        negative_cdr3_series = (
            negative_source.data[negative_source.headers["cdr3"]]
            .sample(n=amount)
            .reset_index(drop=True)
            .rename(negative_source.headers["cdr3"])
        )

        # match with positive epitopes
        negative_df = pd.concat(
            [negative_cdr3_series, self.data[self.headers["epitope_header"]]], axis=1
        )

        # add class labels
        negative_df["y"] = 0

        # merge with positive dataset
        self.data = self.data.append(negative_df).reset_index(drop=True)

        # remove false negatives
        to_do_df = self.data.loc[
            self.data.duplicated(
                subset=[self.headers["cdr3_header"], self.headers["epitope_header"]],
                keep="first",
            )
        ]
        self.data = self.data.drop_duplicates(
            subset=[self.headers["cdr3_header"], self.headers["epitope_header"]],
            keep="first",
        ).reset_index(drop=True)

        # create new negative pairs for any accidental false negatives
        amount = to_do_df.shape[0]
        while amount > 0:
            negative_cdr3_series = (
                negative_source.data[negative_source.headers["cdr3"]]
                .sample(n=amount)
                .reset_index(drop=True)
                .rename("cdr3")
            )

            negative_df = pd.concat(
                [
                    negative_cdr3_series,
                    to_do_df[self.headers["epitope_header"]].reset_index(drop=True),
                ],
                axis=1,
            )

            negative_df["y"] = 0

            self.data = self.data.append(negative_df).reset_index(drop=True)

            to_do_df = self.data.loc[
                self.data.duplicated(
                    subset=[
                        self.headers["cdr3_header"],
                        self.headers["epitope_header"],
                    ],
                    keep="first",
                )
            ]

            self.data = self.data.drop_duplicates(
                subset=[self.headers["cdr3_header"], self.headers["epitope_header"]],
                keep="first",
            ).reset_index(drop=True)

    def generate_negatives(self):
        # generate negative pairs from list of all cdr3s in positive pairs
        shuffled_pairs = [
            sample_pairs(
                cdr3=cdr3,
                df=self.data,
                cdr3_column=self.headers["cdr3_header"],
                epitope_column=self.headers["epitope_header"],
            )
            for cdr3 in self.data[self.headers["cdr3_header"]]
        ]

        # convert list of tuples into dataframe
        shuffled_df = pd.DataFrame(
            shuffled_pairs,
            columns=[self.headers["cdr3_header"], self.headers["epitope_header"]],
        )

        # add negative class labels
        shuffled_df["y"] = 0

        # merge with original positive data, ensuring positives are kept at the top of the dataframe
        self.data = self.data.append(shuffled_df).reset_index(drop=True)

        # extract duplicates
        # NOTE: because the sampling approach ensures that accidental duplicates of positive pairs (i.e. false negatives)
        #       never occur, these will all be accidental duplicate samples of negatives. Therefor, keep="first" is redundant,
        #       but it would result in the positive examples being stored in the dataframe (marking the first (=positives) as True).
        #       This is kept here for historical purposes, because before the epitope was supplied alongside the cdr3 in a zip operation
        #       during sample generation, and the associated positive epitope was used for exclusion purposes.
        to_do_df = self.data.loc[
            self.data.duplicated(
                subset=[self.headers["cdr3_header"], self.headers["epitope_header"]],
                keep="first",
            )
        ]

        # remove duplicates from merged dataframe
        self.data = self.data.drop_duplicates(
            subset=[self.headers["cdr3_header"], self.headers["epitope_header"]],
            keep="first",  # This should not be required, see previous NOTE: always keep the original positive examples when duplicates occur across pos/neg, i.e. remove false negatives
        ).reset_index(drop=True)

        # remove NaN to deal with any possible universal cdr3s
        self.data = self.data.dropna(axis=0, how="any")

        # log difference in positive and negative examples
        # logger = logging.getLogger(__name__)
        # logger.info(
        #     "Negative examples were generated over the entire supplied dataset of positive pairs. The number of positive and negative sequence pairs generated through shuffling and after filtering out duplicates and false negatives is:"
        # )
        # logger.info(self.data["y"].value_counts())

        while to_do_df.shape[0] > 0:
            shuffled_pairs = [
                sample_pairs(
                    cdr3=cdr3,
                    df=self.data,
                    cdr3_column=self.headers["cdr3_header"],
                    epitope_column=self.headers["epitope_header"],
                )
                for cdr3 in to_do_df[self.headers["cdr3_header"]]
            ]
            shuffled_df = pd.DataFrame(
                shuffled_pairs,
                columns=[self.headers["cdr3_header"], self.headers["epitope_header"]],
            )
            shuffled_df["y"] = 0
            self.data = self.data.append(shuffled_df).reset_index(drop=True)
            to_do_df = self.data.loc[
                self.data.duplicated(
                    subset=[
                        self.headers["cdr3_header"],
                        self.headers["epitope_header"],
                    ],
                    keep="first",
                )
            ]
            self.data = self.data.drop_duplicates(
                subset=[self.headers["cdr3_header"], self.headers["epitope_header"]],
                keep="first",
            ).reset_index(drop=True)

            # remove NaN to deal with any possible universal cdr3s
            self.data = self.data.dropna(axis=0, how="any")
