import logging
from typing import Optional, Tuple

from Bio.SubsMat import MatrixInfo
import numpy as np
import pandas as pd
import tensorflow as tf

from src.definitions.amino_acid_properties import AMINO_ACIDS
from src.processing.data_stream import DataStream
from src.processing.negative_sampler import add_negatives
from src.processing.stream import TransformStream
from src.processing.zipper import Zipper


def separated_input_dataset_generator(
    data_stream: DataStream,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
    neg_shuffle: bool = True,
    export_path: Optional[str] = None,
):
    """Create a tensorflow dataset with positive and negative blosum-encoded arrays.

    Can optionally export the positive and generated negative sequence pairs
    to a csv file.

    Parameters
    ----------
    data_stream : DataStream
        A DataStream of positive labeled cdr3-epitope sequence pairs. Expected fromat ( ("CDR3","EPITOPE"), 1)
    cdr3_range : Tuple[int, int]
        The minimum and maximum desired cdr3 sequence length.
    epitope_range : Tuple[int, int]
        The minimum and maximum desired epitope sequence length.
    neg_shuffle : bool
        Whether to create negatives by shuffling/sampling, by default True.
        NOTE: Should always be set to False when evaluating a dataset that already contains negatives.
    export_path : Optional[str], optional
        If supplied, the train/test datasets will be saved to the data/processed directory under this name as a csv file with both positive and negative sequences, by default None.

    Returns
    -------
    tf.data.Dataset
        A tensorflow DataSet, ready to be used as input for a model.
        NOTE: should still be shuffled and batched.
    """
    logger = logging.getLogger(__name__)

    # store maximum cdr3 and epitope sequence length as padding limits
    cdr3_max = cdr3_range[1]
    epitope_max = epitope_range[1]

    # create dataframe for export (and shuffling)
    df = pd.DataFrame(data_stream, columns=["seq", "y"])
    df["cdr3"], df["antigen.epitope"] = zip(*df.seq)
    df = df.drop("seq", axis=1)

    # generate negatives through shuffling if negative reference set was not provided and shuffling did not happen on the entire dataset
    if neg_shuffle:
        df = add_negatives(df)

    # export dataset with sequences and labels as csv
    if export_path:
        logger.info(f"Saving train/test fold in: {export_path}")
        df.to_csv(export_path, sep=";", index=False)

    # Combine sequences and labels into DataStream again to utilise image generation functionality
    zipped = Zipper(
        DataStream(zip(df["cdr3"], df["antigen.epitope"])), DataStream(df["y"])
    )

    # create BLOSUM encoded arrays
    blosum_encoding = BlosumImageGenerator(zipped)

    # pad with zero-columns
    blosum_padding = BlosumPadding(blosum_encoding, cdr3_max, epitope_max, pad_value=0)

    # split stream back into separate sequence and label tuples for export to tf DataSet
    x, y = zip(*blosum_padding)
    x_cdr3, x_epitope = zip(*x)

    # convert into tf DataSet
    dataset = tf.data.Dataset.from_tensor_slices(
        (np.array(x_cdr3), np.array(x_epitope), np.array(y))
    )

    return dataset


class BlosumImageGenerator(TransformStream):
    def __init__(self, stream):
        super().__init__(stream)
        self.features = dict()
        for aa in AMINO_ACIDS:
            self.features[aa] = [self._get_matrix_entry(aa, x) for x in AMINO_ACIDS]

    def _get_matrix_entry(self, aa1, aa2):
        i = MatrixInfo.blosum50.get((aa1, aa2))
        if i is not None:
            return i
        else:
            return MatrixInfo.blosum50.get((aa2, aa1))

    def transform(self, item, *args, **kwargs):
        x, y = item
        X = tuple(self.convert(xx) for xx in x)
        return X, y

    def convert(self, sequence):
        array = np.array([self.features[aa] for aa in sequence])
        return array


class BlosumPadding(TransformStream):
    def __init__(self, stream, max_length_cdr3, max_length_epitope, pad_value=0):
        super().__init__(stream)
        self.max_length_cdr3 = max_length_cdr3
        self.max_length_epitope = max_length_epitope
        self.pad_value = pad_value

    def transform(self, item, *args, **kwargs):
        sequence_tuple, label = item
        cdr3_array, epitope_array = sequence_tuple

        cdr3_padding = ((0, self.max_length_cdr3 - cdr3_array.shape[0]), (0, 0))
        epitope_padding = (
            (0, self.max_length_epitope - epitope_array.shape[0]),
            (0, 0),
        )

        padded_cdr3_array = np.pad(
            cdr3_array, cdr3_padding, mode="constant", constant_values=self.pad_value
        )
        padded_epitope_array = np.pad(
            epitope_array,
            epitope_padding,
            mode="constant",
            constant_values=self.pad_value,
        )

        return (padded_cdr3_array, padded_epitope_array), label