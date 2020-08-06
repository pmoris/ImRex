import logging
from pathlib import Path
from typing import Optional, Tuple

from Bio.SubsMat import MatrixInfo
import numpy as np
import pandas as pd
import tensorflow as tf

from src.definitions.amino_acid_properties import AMINO_ACIDS
from src.processing.data_stream import DataStream
from src.processing.negative_sampler import add_negatives, augment_negatives
from src.processing.stream import TransformStream
from src.processing.zipper import Zipper


def separated_input_dataset_generator(
    data_stream: DataStream,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
    neg_shuffle: bool = True,
    full_dataset_path: Optional[Path] = None,
    epitope_ratio: bool = False,
    export_path: Optional[str] = None,
    neg_augment: Optional[str] = None,
    augment_amount: Optional[int] = None,
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
    full_dataset_path : Path
        The entire cdr3-epitope dataset, before splitting into folds, restricting length or downsampling. Used to avoid
        generating false negatives. Should only contain positive values.
    epitope_ratio : boolean
        When false, samples an epitope for each CDR3 sequence in the
        proportionally to its occurrence in the other epitope pairs. Does not
        preserve the ratio of positives and negatives within each epitope,
        but does result in every CDR3 sequence having exactly 1 positive and negative.
        When true, samples a set of CDR3 sequences with from the unique list of CDR3s
        for each epitope observation (per epitope), i.e. preserves exact ratio of positives and
        negatives for each epitope, at the expense of some CDR3s appearing more than once
        among the negatives and others only in positives pairs.
    export_path: Optional[str], optional
        If supplied, the train/test datasets will be saved to the data/processed directory under this name as a csv file with both positive and negative sequences, by default None.
    neg_augment: Optional[str], optional
        If supplied, provided the filepath to a negative reference set of cdr3 sequences, used for augmenting additional negatives, by default None.
    augment_amount: Optional[int], optional
        The amount of negatives to augment.

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
        assert (
            full_dataset_path
        ), "The path to the full dataset should be supplied when generating negatives through shuffling."
        df = add_negatives(
            df=df, full_dataset_path=full_dataset_path, epitope_ratio=epitope_ratio
        )

    # optionally augment with additional negative reference pairs
    if neg_augment and augment_amount:
        logger.info(
            f"Augmenting {augment_amount} negatives from negative reference set {Path(neg_augment).absolute()}."
        )
        df = augment_negatives(
            negative_source=neg_augment,
            df=df,
            cdr3_range=cdr3_range,
            amount=augment_amount,
        )

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
        ((np.array(x_cdr3), np.array(x_epitope)), np.array(y))
    )

    return dataset


class BlosumImageGenerator(TransformStream):
    def __init__(self, stream, has_label=True):
        super().__init__(stream)
        self.features = dict()
        for aa in AMINO_ACIDS:
            self.features[aa] = [self._get_matrix_entry(aa, x) for x in AMINO_ACIDS]
        self.has_label = has_label

    def _get_matrix_entry(self, aa1, aa2):
        i = MatrixInfo.blosum50.get((aa1, aa2))
        if i is not None:
            return i
        else:
            return MatrixInfo.blosum50.get((aa2, aa1))

    def transform(self, item, *args, **kwargs):
        if self.has_label:
            x, y = item
            X = tuple(self.convert(xx) for xx in x)
            return X, y
        else:
            X = tuple(self.convert(xx) for xx in item)
            return X

    def convert(self, sequence):
        array = np.array([self.features[aa] for aa in sequence])
        return array


class BlosumPadding(TransformStream):
    def __init__(
        self, stream, max_length_cdr3, max_length_epitope, pad_value=0, has_label=True
    ):
        super().__init__(stream)
        self.max_length_cdr3 = max_length_cdr3
        self.max_length_epitope = max_length_epitope
        self.pad_value = pad_value
        self.has_label = has_label

    def transform(self, item, *args, **kwargs):
        if self.has_label:
            sequence_tuple, label = item
        else:
            sequence_tuple = item

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

        if self.has_label:
            return (padded_cdr3_array, padded_epitope_array), label
        else:
            return (padded_cdr3_array, padded_epitope_array)

