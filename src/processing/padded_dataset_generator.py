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
    batch_size: int,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
    inverse_map: Optional[InverseMap] = NoOp(),
    negative_stream: Optional[DataStream] = None,
):
    logger = logging.getLogger(__name__)

    width = cdr3_range[1]
    height = epitope_range[1]

    # filter on sequences on size
    size_filter = SizeFilter(data_stream, cdr3_range, epitope_range)
    logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
    logger.info(f"Filtered epitope sequences to length: {epitope_range}")

    # split into 3 streams: positive, negative and positive filter
    input1, input2, input3 = tee(size_filter, amount=3)

    input1 = inverse_map.input(input1)
    pos_images = ImageGenerator(input1, feature_builder)
    pos_out = ImagePadding(pos_images, width, height, pad_value=0)
    pos_out = inverse_map.output(pos_out)

    # transform into np.array
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

    # sample random cdr3 and epitopes once and retain this throughout epochs
    cdr3, epitope = np.array(list(cdr3)), np.array(list(epitope))
    np.random.shuffle(cdr3)
    np.random.shuffle(epitope)
    cdr3, epitope = DataStream(cdr3), DataStream(epitope)
    zipper = Zipper(cdr3, epitope)

    # remove matches with existing positive pairs
    pos_filter = PositiveFilter(zipper, positive_items=input3, has_label=False)

    # add negative class labels
    negative_labeler = Labeler(pos_filter, 0)

    negative_labeler = inverse_map.input(negative_labeler)
    neg_image_gen = ImageGenerator(negative_labeler, feature_builder)
    neg_padding = ImagePadding(neg_image_gen, width, height, pad_value=0)
    neg_padding = inverse_map.output(neg_padding)

    x_neg, y_neg = zip(*neg_padding)
    logger.info(f"Using {len(x_neg)} negative examples.")

    x, y = np.array(x_pos + x_neg), np.array(y_pos + y_neg)
    dataset = tf.data.Dataset.from_tensor_slices((x, y))
    dataset = dataset.shuffle(
        buffer_size=len(x), seed=42, reshuffle_each_iteration=True
    ).batch(batch_size)

    return dataset
