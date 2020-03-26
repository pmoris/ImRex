import logging
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import tensorflow as tf

from src.bio.feature_builder import FeatureBuilder
from src.processing.data_stream import DataStream
from src.processing.image_generator import ImageGenerator
from src.processing.image_padding import ImagePadding
from src.processing.inverse_map import InverseMap, NoOp
from src.processing.negative_sampler import add_negatives
from src.processing.zipper import Zipper


def padded_dataset_generator(
    data_stream: DataStream,
    feature_builder: FeatureBuilder,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
    inverse_map: Optional[InverseMap] = NoOp(),
    neg_shuffle: bool = True,
    export_path: Optional[str] = None,
) -> tf.data.Dataset:
    """Create a tensorflow dataset with positive and negative 2d interaction map arrays.

    Can optionally export the positive and generated negative sequence pairs
    to a csv file.

    Parameters
    ----------
    data_stream : DataStream
        A DataStream of positive labeled cdr3-epitope sequence pairs.
    feature_builder : FeatureBuilder
        A FeatureBuilder object that can convert the sequences into pairwise interaction arrays.
    cdr3_range : Tuple[int, int]
        The minimum and maximum desired cdr3 sequence length.
    epitope_range : Tuple[int, int].
        The minimum and maximum desired epitope sequence length.
    inverse_map : Optional[InverseMap], optional.
        An inverse map for retrieving the sequences associated with an array, by default NoOp().
    neg_shuffle : bool
        Whether to create negatives by shuffling/sampling, by default True.
    export_path: Optional[str], optional.
        If supplied, the train/test datasets will be saved to the data/processed directory under this name as a csv file with both positive and negative sequences, by default None.

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
        df = add_negatives(df)

    # export dataset with sequences and labels as csv
    if export_path:
        logger.info(f"Saving train/test fold in: {export_path}")
        df.to_csv(export_path, sep=";", index=False)

    # Combine sequences and labels into DataStream again to utilise image generation functionality
    zipped = Zipper(
        DataStream(zip(df["cdr3"], df["antigen.epitope"])), DataStream(df["y"])
    )

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
