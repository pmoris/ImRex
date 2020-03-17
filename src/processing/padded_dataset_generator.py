import logging
from typing import Optional, Tuple

import numpy as np
import tensorflow as tf

from src.bio.feature_builder import FeatureBuilder
from src.processing.data_stream import DataStream
from src.processing.filter import PositiveFilter, SizeFilter
from src.processing.image_generator import ImageGenerator
from src.processing.image_padding import ImagePadding
from src.processing.inverse_map import InverseMap, NoOp
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.tee import tee
from src.processing.zipper import unzipper, Zipper


def padded_dataset_generator(
    data_stream: DataStream,
    feature_builder: FeatureBuilder,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
    inverse_map: Optional[InverseMap] = NoOp(),
    negative_stream: Optional[DataStream] = None,
) -> tf.data.Dataset:
    """Create a tensorflow dataset with positive and negative 2d arrays.

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
    negative_stream : Optional[DataStream], optional
        An optional stream of reference cdr3 sequences, without epitope or label, by default None.

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

    # filter sequences on size
    size_filter = SizeFilter(data_stream, cdr3_range, epitope_range)
    logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
    logger.info(f"Filtered epitope sequences to length: {epitope_range}")

    # replicate FilterStream into 3 DataStreams: for positive examples, negative examples and positive filter
    # NOTE: SizeFilter is exhausted upon iteration/instantiation, so it is copied back into a DataStream via tee.
    input1, input2, input3 = tee(size_filter, amount=3)

    # create positive 2d input arrays and allow reverse lookup back to sequences
    input1 = inverse_map.input(input1)
    pos_images = ImageGenerator(input1, feature_builder)
    pos_out = ImagePadding(pos_images, width, height, pad_value=0)
    pos_out = inverse_map.output(pos_out)

    # transform into np.array for later conversion into tf DataSet
    x_pos, y_pos = zip(*pos_out)
    logger.info(f"Using {len(x_pos)} positive examples.")

    # Create negative stream

    # remove class labels from stream
    label_trimmer = LabelTrimmer(input2)
    # split cdr3 and epitope sequences into separate DataStreams
    cdr3, epitope = unzipper(label_trimmer)
    # if reference cdr3 was provided use it
    if negative_stream:
        cdr3 = SizeFilter(negative_stream, cdr3_range, has_label=False)

    # Create permuted cdr3 and epitopes sequence pairs (retained throughout consecutive iterations of DataSet, and thus epochs).

    # convert DataStreams to arrays
    cdr3, epitope = np.array(list(cdr3)), np.array(list(epitope))

    # create permutations
    cdr3_permuted = np.random.permutation(cdr3)
    epitope_permuted = np.random.permutation(epitope)

    # convert back into DataStreams
    cdr3_permuted_stream, epitope_permuted_stream = (
        DataStream(cdr3_permuted),
        DataStream(epitope_permuted),
    )

    # combine sequences
    zipper = Zipper(cdr3_permuted_stream, epitope_permuted_stream)
    # remove matches with existing positive pairs
    pos_filter = DataStream(
        PositiveFilter(zipper, positive_items=input3, has_label=False)
    )  # NOTE: cast into DataStream so that its length becomes available

    # create additional random pairs to even out the number of positive and negative examples
    while len(pos_filter) < len(x_pos):
        # calculate how many new pairs are still required
        amount = len(x_pos) - len(pos_filter)
        pos_filter = augment_pairs(
            data_stream=pos_filter,
            amount=amount,
            positive_filter_set=input3,
            negative_stream=negative_stream,
            cdr3_range=cdr3_range,
        )

    # add negative class labels
    negative_labeler = Labeler(pos_filter, 0)

    # create negative 2d input arrays and allow reverse lookup back to sequences
    negative_labeler = inverse_map.input(negative_labeler)
    neg_image_gen = ImageGenerator(negative_labeler, feature_builder)
    neg_padding = ImagePadding(neg_image_gen, width, height, pad_value=0)
    neg_padding = inverse_map.output(neg_padding)

    x_neg, y_neg = zip(*neg_padding)
    logger.info(f"Using {len(x_neg)} negative examples.")

    x, y = np.array(x_pos + x_neg), np.array(y_pos + y_neg)
    dataset = tf.data.Dataset.from_tensor_slices((x, y))

    return dataset


def augment_pairs(
    data_stream: DataStream,
    amount: int,
    positive_filter_set: DataStream,
    cdr3_range: Tuple[int, int] = None,
    negative_stream: Optional[DataStream] = None,
):
    # split DataStream into two Streams
    cdr3_permuted_stream, epitope_permuted_stream = unzipper(data_stream)

    # if reference cdr3 was provided use it
    if negative_stream:
        # assign negative_stream to a new DataStream object to reset exhausted iterator
        negative_stream = DataStream(negative_stream)
        cdr3_permuted_stream = SizeFilter(negative_stream, cdr3_range, has_label=False)

    # convert to np.arrays to allow concatenation
    cdr3_permuted, epitope_permuted = (
        np.array(list(cdr3_permuted_stream)),
        np.array(list(epitope_permuted_stream)),
    )

    # create new permutations (it does not matter if the already permuted arrays or the original arrays are shuffled again)
    cdr3_permuted_new = np.random.permutation(cdr3_permuted)[:amount]
    epitope_permuted_new = np.random.permutation(epitope_permuted)[:amount]

    # concatenate old and new arrays and convert to DataStreams
    cdr3_permuted_stream = DataStream(
        np.concatenate([cdr3_permuted, cdr3_permuted_new])
    )
    epitope_permuted_stream = DataStream(
        np.concatenate([epitope_permuted, epitope_permuted_new])
    )

    # combine cdr3 and epitope sequences
    zipper = Zipper(cdr3_permuted_stream, epitope_permuted_stream)

    # remove matches with existing positive pairs
    pos_filter = DataStream(
        PositiveFilter(zipper, positive_items=positive_filter_set, has_label=False)
    )

    return pos_filter
