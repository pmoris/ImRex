import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import tensorflow as tf

from src.bio.feature_builder import FeatureBuilder
from src.processing.data_stream import DataStream
from src.processing.image_generator import ImageGenerator
from src.processing.image_padding import ImagePadding
from src.processing.inverse_map import InverseMap, NoOp
from src.processing.negative_sampler import add_negatives, augment_negatives
from src.processing.zipper import Zipper


def padded_dataset_generator(
    data_stream: DataStream,
    feature_builder: FeatureBuilder,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
    inverse_map: Optional[InverseMap] = NoOp(),
    neg_shuffle: bool = True,
    full_dataset_path: Optional[Path] = None,
    epitope_ratio: bool = False,
    export_path: Optional[str] = None,
    neg_augment: Optional[str] = None,
    augment_amount: Optional[int] = None,
    reverse_augment: Optional[bool] = False,
) -> tf.data.Dataset:
    """Create a tensorflow dataset with positive and negative 2d interaction map arrays.

    Can optionally export the positive and generated negative sequence pairs
    to a csv file.

    Parameters
    ----------
    data_stream : DataStream
        A DataStream of positive labeled cdr3-epitope sequence pairs. Expected fromat ( ("CDR3","EPITOPE"), 1)
    feature_builder : FeatureBuilder
        A FeatureBuilder object that can convert the sequences into pairwise interaction arrays.
    cdr3_range : Tuple[int, int]
        The minimum and maximum desired cdr3 sequence length.
    epitope_range : Tuple[int, int]
        The minimum and maximum desired epitope sequence length.
    inverse_map : Optional[InverseMap], optional
        An inverse map for retrieving the sequences associated with an array, by default NoOp().
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

    # store maximum cdr3 and epitope sequence length as width and height of the 2d arrays
    width = cdr3_range[1]
    height = epitope_range[1]

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

    # augment data by reversing sequences
    if reverse_augment:
        df_rev_cdr3 = df.copy()
        df_rev_cdr3["cdr3"] = df_rev_cdr3["cdr3"].map(lambda x: x[::-1])

        df_rev_epitope = df.copy()
        df_rev_epitope["antigen.epitope"] = df_rev_epitope["antigen.epitope"].map(
            lambda x: x[::-1]
        )

        df_rev_both = df.copy()
        df_rev_both["cdr3"] = df_rev_both["cdr3"].map(lambda x: x[::-1])
        df_rev_both["antigen.epitope"] = df_rev_both["antigen.epitope"].map(
            lambda x: x[::-1]
        )

        df = pd.concat([df, df_rev_cdr3, df_rev_epitope, df_rev_both])
        df = df.reset_index(drop=True)

    # export dataset with sequences and labels as csv
    if export_path:
        logger.info(f"Saving train/test fold in: {export_path}")
        df.to_csv(export_path, sep=";", index=False)

    # Combine sequences and labels into DataStream again to utilise image generation functionality
    zipped = Zipper(
        DataStream(zip(df["cdr3"], df["antigen.epitope"])), DataStream(df["y"])
    )

    # # augment data by reversing sequences
    # # does not work because zipped is exhausted on the first pass
    # if reverse_augment:
    #     from itertools import chain

    #     rev_cdr3 = Reverser(zipped, target="cdr3")
    #     rev_epitope = Reverser(zipped, target="epitope")
    #     rev_both = Reverser(zipped, target=None)
    #     zipped = DataStream(chain(rev_cdr3, rev_epitope, rev_both))
    #     logger.info(
    #         f"Augmenting data by reversing sequence direction. Total number of training samples is {len(zipped)}."
    #     )

    # create 2d input arrays and allow reverse lookup back to sequences
    zipped = inverse_map.input(zipped)
    image_gen = ImageGenerator(zipped, feature_builder)
    image_padding = ImagePadding(image_gen, width, height, pad_value=0)
    image_stream = inverse_map.output(image_padding)

    # split stream back into separate sequence and label tuples for export to tf DataSet
    x, y = zip(*image_stream)

    # convert into tf DataSet
    dataset = tf.data.Dataset.from_tensor_slices((np.array(x), np.array(y)))

    return dataset
