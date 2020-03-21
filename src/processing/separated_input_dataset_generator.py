import logging
from typing import Optional, Tuple

from Bio.SubsMat import MatrixInfo
import numpy as np
import pandas as pd
import tensorflow as tf

from src.definitions.amino_acid_properties import AMINO_ACIDS
from src.processing.data_stream import DataStream
from src.processing.dataset_export import export_data
from src.processing.filter import PositiveFilter, SizeFilter
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.stream import TransformStream
from src.processing.zipper import unzipper, Zipper


def separated_input_dataset_generator(
    data_stream: DataStream,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
    negative_ref_stream: Optional[DataStream] = None,
    export_path: Optional[str] = None,
):
    """Create a tensorflow dataset with positive and negative blosum-encoded arrays.

    Can optionally export the positive and generated negative sequence pairs
    to a csv file.

    Parameters
    ----------
    data_stream : DataStream
        A DataStream of positive labeled cdr3-epitope sequence pairs.
    cdr3_range : Tuple[int, int]
        The minimum and maximum desired cdr3 sequence length.
    epitope_range : Tuple[int, int]
        The minimum and maximum desired epitope sequence length.
    negative_ref_stream : Optional[DataStream], optional
        An optional stream of reference cdr3 sequences, without epitope or label, by default None.
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

    # create positive sequences

    # filter sequences on size
    size_filter = SizeFilter(data_stream, cdr3_range, epitope_range)
    logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
    logger.info(f"Filtered epitope sequences to length: {epitope_range}")

    # NOTE: SizeFilter is exhausted upon iteration/instantiation, so it is copied back into a DataStream to be re-used for negative example generation
    positive_stream = DataStream(size_filter)

    # split positive sequences and labels into separate tuples for later combining and export
    x_pos_seq, y_pos_seq = zip(*positive_stream)

    # Create negatives by permuting cdr3 and epitopes sequence pairs
    # These are retained throughout consecutive iterations of the tf DataSet, and thus epochs.

    # remove class labels
    label_trimmer = LabelTrimmer(positive_stream)

    # split DataStream into two Streams
    cdr3_stream, epitope_stream = unzipper(label_trimmer)

    # # if reference cdr3 was provided use it
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
    # logger.info(f"Using {len(x_neg_seq)} negative examples.")

    # combine positive and negative examples into arrays with sequences and labels
    x_seq, y_seq = np.array(x_pos_seq + x_neg_seq), np.array(y_pos_seq + y_neg_seq)

    # export dataset with sequences and labels as csv
    if export_path:
        export_data(x_seq, y_seq, export_path)

    # Combine sequences and labels into DataStream again to utilise image generation functionality
    zipped = Zipper(DataStream(x_seq), DataStream(y_seq))

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
