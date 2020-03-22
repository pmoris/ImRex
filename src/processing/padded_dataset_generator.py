import logging
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import tensorflow as tf

from src.bio.feature_builder import FeatureBuilder
from src.processing.data_stream import DataStream
from src.processing.dataset_export import export_data
from src.processing.filter import PositiveFilter, SizeFilter
from src.processing.image_generator import ImageGenerator
from src.processing.image_padding import ImagePadding
from src.processing.inverse_map import InverseMap, NoOp
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.zipper import unzipper, Zipper


def padded_dataset_generator(
    data_stream: DataStream,
    feature_builder: FeatureBuilder,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
    inverse_map: Optional[InverseMap] = NoOp(),
    negative_ref_stream: Optional[DataStream] = None,
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
    negative_ref_stream : Optional[DataStream], optional
        An optional stream of reference cdr3 sequences, without epitope or label, by default None.
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

    # create positive sequences

    # filter sequences on size
    size_filter = SizeFilter(data_stream, cdr3_range, epitope_range)
    logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
    logger.info(f"Filtered epitope sequences to length: {epitope_range}")

    # NOTE: SizeFilter is exhausted upon iteration/instantiation, so it is copied back into a DataStream to be re-used for negative example generation
    positive_stream = DataStream(size_filter)

    # split positive sequences and labels into separate tuples for later combining and export
    x_pos_seq, y_pos_seq = zip(*positive_stream)
    logger.info(f"Using {len(x_pos_seq)} positive examples.")

    # Create negatives by permuting cdr3 and epitopes sequence pairs
    # These are retained throughout consecutive iterations of the tf DataSet, and thus epochs.

    # remove class labels
    label_trimmer = LabelTrimmer(positive_stream)

    # split DataStream into two Streams
    cdr3_stream, epitope_stream = unzipper(label_trimmer)

    # if reference cdr3 was provided use it
    if negative_ref_stream:
        cdr3_stream = SizeFilter(negative_ref_stream, cdr3_range, has_label=False)

    # convert DataStreams to np.arrays to allow shuffling and concatenation
    cdr3_array, epitope_array = (
        np.array(list(cdr3_stream)),
        np.array(list(epitope_stream)),
    )

    # create permutations
    cdr3_permuted = np.random.permutation(cdr3_array)
    epitope_permuted = np.random.permutation(epitope_array)

    # convert back into DataStreams to utilise zipper/filter functionality
    cdr3_permuted_stream, epitope_permuted_stream = (
        DataStream(cdr3_permuted),
        DataStream(epitope_permuted),
    )

    # combine sequences
    zipper = Zipper(cdr3_permuted_stream, epitope_permuted_stream)

    # remove matches with existing positive pairs
    pos_filter = DataStream(
        PositiveFilter(zipper, positive_items=positive_stream, has_label=False)
    )  # NOTE: cast into DataStream so that its length becomes available

    # create additional random pairs to even out the number of positive and negative examples
    while len(pos_filter) < len(positive_stream):

        # calculate how many new pairs are still required
        amount = len(positive_stream) - len(pos_filter)

        pos_filter = augment_pairs(
            permuted_stream=pos_filter,
            amount=amount,
            original_cdr3_array=cdr3_array,
            original_epitope_array=epitope_array,
            positive_filter_set=positive_stream,
            cdr3_range=cdr3_range,
        )

    # add negative class labels
    negative_labeler = Labeler(pos_filter, 0)

    # split negative sequences and labels into separate tuples for later combining and export
    x_neg_seq, y_neg_seq = zip(*negative_labeler)
    logger.info(f"Using {len(x_neg_seq)} negative examples.")

    # combine positive and negative examples into arrays with sequences and labels
    x_seq, y_seq = np.array(x_pos_seq + x_neg_seq), np.array(y_pos_seq + y_neg_seq)

    # export dataset with sequences and labels as csv
    if export_path:
        export_data(x_seq, y_seq, export_path)

    # Combine sequences and labels into DataStream again to utilise image generation functionality
    zipped = Zipper(DataStream(x_seq), DataStream(y_seq))

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


def augment_pairs(
    permuted_stream: DataStream,
    amount: int,
    original_cdr3_array: np.array,
    original_epitope_array: np.array,
    positive_filter_set: DataStream,
    cdr3_range: Tuple[int, int] = None,
) -> DataStream:
    # split permuted DataStream into two Streams
    cdr3_stream, epitope_stream = unzipper(permuted_stream)

    # convert DataStreams to np.arrays to allow concatenation
    cdr3_array, epitope_array = (
        np.array(list(cdr3_stream)),
        np.array(list(epitope_stream)),
    )

    # create new permutations from original arrays
    cdr3_permuted_new = np.random.permutation(original_cdr3_array)[:amount]
    epitope_permuted_new = np.random.permutation(original_epitope_array)[:amount]

    # concatenate old and new arrays and convert to DataStreams
    cdr3_permuted_stream = DataStream(np.concatenate([cdr3_array, cdr3_permuted_new]))
    epitope_permuted_stream = DataStream(
        np.concatenate([epitope_array, epitope_permuted_new])
    )

    # combine cdr3 and epitope sequences
    zipper = Zipper(cdr3_permuted_stream, epitope_permuted_stream)

    # remove duplicates due to repeated shuffling
    deduplicated_df = pd.DataFrame(zipper).drop_duplicates()
    deduplicated_stream = DataStream(deduplicated_df.to_numpy().tolist())

    # remove matches with existing positive pairs
    pos_filter = DataStream(
        PositiveFilter(
            deduplicated_stream, positive_items=positive_filter_set, has_label=False
        )
    )

    return pos_filter
