import numpy as np

from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.control_cdr3_source import ControlCDR3Source
from src.data.vdjdb_source import VdjdbSource
from src.processing.data_stream import DataStream
from src.processing.filter import PositiveFilter
from src.processing.labeler import LabelTrimmer
from src.processing.padded_batch_generator import padded_batch_generator
from src.processing.padded_dataset_generator import (
    augment_pairs,
    padded_dataset_generator,
)
from src.processing.zipper import unzipper, Zipper


features_list = parse_features("hydrophob,isoelectric,mass,hydrophil,charge")
operator = parse_operator("absdiff")
feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)


def test_tf_dataset_shuffle():
    """Test whether sequential iteration through tf.data.DataSet objects results in different batches.

    E.g. during every epoch, the generated batches should ideally be different every time.

    When using a DataSet object that was created from a NumPy array, this behaviour is controlled by the
    shuffle(buffer_size=len(x), seed=42, reshuffle_each_iteration=True)
    reshuffle_each_iteration argument.

    However, when using a custom generator, the behaviour is controlled by the underlying generator code.
    """
    data_stream = DataStream(
        VdjdbSource(
            filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
            headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
        )
    )

    tf_dataset = padded_batch_generator(
        data_stream=data_stream,
        feature_builder=feature_builder,
        neg_ratio=0.5,
        cdr3_range=(10, 20),
        epitope_range=(8, 11),
        negative_ref_stream=None,
    )

    tf_dataset = tf_dataset.batch(10)

    iteration_1 = [i[1] for i in list(tf_dataset.as_numpy_iterator())]
    iteration_2 = [i[1] for i in list(tf_dataset.as_numpy_iterator())]

    # list of arrays where each array is a batch. equality of np.arrays results in array of bools per comparison
    # [array([False, False, False, False,  True]), array([False, False, False, False, False]), ...
    # print([i == j for i, j in zip(iteration_1, iteration_2)])

    # [False, False, False, ...] per array/batch
    # print([all(a) for a in [i == j for i, j in zip(iteration_1, iteration_2)]])

    # as long as not all arrays/batches are completely identical, there is shuffling happening between iterations
    assert not all([all(a) for a in [i == j for i, j in zip(iteration_1, iteration_2)]])


def test_tf_dataset_shuffle_array():
    """Same as above, but create the tf.data.DataSet object from numpy arrays. """
    # NOTE: DataStream needs to be re-created, because it will be exhausted by previous test otherwise.
    data_stream = DataStream(
        VdjdbSource(
            filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
            headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
        )
    )

    tf_dataset = padded_dataset_generator(
        data_stream=data_stream,
        feature_builder=feature_builder,
        cdr3_range=(10, 20),
        epitope_range=(8, 11),
        negative_ref_stream=None,
    )

    tf_dataset = tf_dataset.shuffle(
        buffer_size=len(data_stream), seed=42, reshuffle_each_iteration=True
    ).batch(10)

    iteration_1 = [i[1] for i in list(tf_dataset.as_numpy_iterator())]
    iteration_2 = [i[1] for i in list(tf_dataset.as_numpy_iterator())]

    assert not all([all(a) for a in [i == j for i, j in zip(iteration_1, iteration_2)]])


def test_tf_dataset_shuffle_array_neg_ref():
    """Same as above, but create the tf.data.DataSet object from numpy arrays. """
    # NOTE: DataStream needs to be re-created, because it will be exhausted by previous test otherwise.
    data_stream = DataStream(
        VdjdbSource(
            filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
            headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
        )
    )

    negative_source = ControlCDR3Source(
        filepath=PROJECT_ROOT / "src/tests/test_CDR3_control.tsv",
        min_length=10,
        max_length=20,
    )
    negative_ref_stream = DataStream(negative_source)

    tf_dataset = padded_dataset_generator(
        data_stream=data_stream,
        feature_builder=feature_builder,
        cdr3_range=(10, 20),
        epitope_range=(8, 11),
        negative_ref_stream=negative_ref_stream,
    )

    tf_dataset = tf_dataset.shuffle(
        buffer_size=len(data_stream), seed=42, reshuffle_each_iteration=True
    ).batch(10)

    iteration_1 = [i[1] for i in list(tf_dataset.as_numpy_iterator())]
    iteration_2 = [i[1] for i in list(tf_dataset.as_numpy_iterator())]

    assert not all([all(a) for a in [i == j for i, j in zip(iteration_1, iteration_2)]])


def test_augment_pairs_dedup():
    """
    Check whether duplicates and positives are correctly filtered out.

    The maximum number of unique permutations of 6 pairs is 36.
    Removing the original 6 pairs leaves 30.
    """
    stream = DataStream(
        [
            (("AAAAAAAAAA", "AAAAAAAAAA"), 1),
            (("BBBBBBBBBB", "BBBBBBBBBB"), 1),
            (("C", "C"), 1),
            (("DDDDDDDDDD", "DDDDDDDDDD"), 1),
            (("EEE", "EEEE"), 1),
            (("FFFFFFFFFFFFF", "FFFFFFFFF"), 1),
        ]
    )
    negative_ref_stream = None
    label_trimmer = LabelTrimmer(stream)
    stream = DataStream(label_trimmer)

    cdr3_stream, epitope_stream = unzipper(stream)
    cdr3_array, epitope_array = (
        np.array(list(cdr3_stream)),
        np.array(list(epitope_stream)),
    )

    for _ in range(100):
        stream = augment_pairs(
            permuted_stream=stream,
            amount=len(stream),
            original_cdr3_array=cdr3_array,
            original_epitope_array=epitope_array,
            positive_filter_set=DataStream([]),
            cdr3_range=(0, 30),
            negative_ref_stream=negative_ref_stream,
        )

    assert len(stream) == 36


def test_augment_pairs_dedup_pos():
    """
    Check whether duplicates are correctly filtered out.

    The maximum number of unique permutations of 6 pairs is 36.
    """
    positive_stream = DataStream(
        [
            (("AAAAAAAAAA", "AAAAAAAAAA"), 1),
            (("BBBBBBBBBB", "BBBBBBBBBB"), 1),
            (("C", "C"), 1),
            (("DDDDDDDDDD", "DDDDDDDDDD"), 1),
            (("EEE", "EEEE"), 1),
            (("FFFFFFFFFFFFF", "FFFFFFFFF"), 1),
        ]
    )
    negative_ref_stream = None
    label_trimmer = LabelTrimmer(positive_stream)
    trimmer = DataStream(label_trimmer)

    cdr3_stream, epitope_stream = unzipper(trimmer)

    cdr3_array, epitope_array = (
        np.array(list(cdr3_stream)),
        np.array(list(epitope_stream)),
    )

    cdr3_permuted = np.random.permutation(cdr3_array)
    epitope_permuted = np.random.permutation(epitope_array)

    cdr3_permuted_stream, epitope_permuted_stream = (
        DataStream(cdr3_permuted),
        DataStream(epitope_permuted),
    )

    zipper = Zipper(cdr3_permuted_stream, epitope_permuted_stream)

    pos_filter = DataStream(
        PositiveFilter(zipper, positive_items=positive_stream, has_label=False)
    )

    for _ in range(100):
        pos_filter = augment_pairs(
            permuted_stream=pos_filter,
            amount=len(positive_stream),
            original_cdr3_array=cdr3_array,
            original_epitope_array=epitope_array,
            positive_filter_set=positive_stream,
            cdr3_range=(0, 30),
            negative_ref_stream=negative_ref_stream,
        )

    assert len(pos_filter) == 30
