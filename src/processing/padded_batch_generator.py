import logging
from typing import Optional, Tuple

import tensorflow as tf

from src.bio.feature_builder import FeatureBuilder
from src.processing.data_stream import DataStream
from src.processing.filter import PositiveFilter, SizeFilter
from src.processing.image_generator import ImageGenerator
from src.processing.image_padding import ImagePadding
from src.processing.inverse_map import InverseMap, NoOp
from src.processing.joiner import Joiner
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.sampler import Sampler
from src.processing.swapper import Swapper
from src.processing.tee import tee
from src.processing.zipper import unzipper, Zipper


def padded_batch_generator(
    data_stream: DataStream,
    feature_builder: FeatureBuilder,
    neg_ratio: float,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
    inverse_map: Optional[InverseMap] = NoOp(),
    negative_ref_stream: Optional[DataStream] = None,
    cache_images: bool = True,
    swap: bool = False,
):
    """ Standard PaddedBatchGenerator. """
    logger = logging.getLogger(__name__)

    width = cdr3_range[1]
    height = epitope_range[1]

    # filter on sequences on size
    size_filter = SizeFilter(data_stream, cdr3_range, epitope_range)
    logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
    logger.info(f"Filtered epitope sequences to length: {epitope_range}")

    # split into 3 streams: positive, negative and positive filter
    input1, input2, input3 = tee(size_filter, amount=3)

    # create positive stream
    if not cache_images:
        input1 = Sampler(input1)

    if swap:
        input1 = Swapper(input1)

    input1 = inverse_map.input(input1)
    pos_images = ImageGenerator(input1, feature_builder)
    pos_out = ImagePadding(pos_images, width, height, pad_value=0)
    pos_out = inverse_map.output(pos_out)

    if cache_images:
        pos_out = Sampler(pos_out)

    # Create negative stream

    # remove class labels from stream
    label_trimmer = LabelTrimmer(input2)
    # split cdr3 and epitope sequences into separate DataStreams
    cdr3, epitope = unzipper(label_trimmer)

    # if reference cdr3 was provided use it
    if negative_ref_stream:
        cdr3 = SizeFilter(negative_ref_stream, cdr3_range, has_label=False)

    # sample random cdr3 and epitopes: repeated for every epoch
    sampler1 = Sampler(cdr3, infinite=True)
    sampler2 = Sampler(epitope, infinite=True)
    zipper = Zipper(sampler1, sampler2)

    # remove matches with existing positive pairs
    pos_filter = PositiveFilter(zipper, positive_items=input3, has_label=False)

    # add negative class labels
    negative_labeler = Labeler(pos_filter, 0)

    if swap:
        negative_labeler = Swapper(negative_labeler)

    negative_labeler = inverse_map.input(negative_labeler)
    neg_image_gen = ImageGenerator(negative_labeler, feature_builder)
    neg_padding = ImagePadding(neg_image_gen, width, height, pad_value=0)
    neg_padding = inverse_map.output(neg_padding)

    # join positive and negatives
    joiner = Joiner(pos_out, neg_padding, 1 - neg_ratio)

    dataset = tf.data.Dataset.from_generator(
        generator=get_generator(joiner, len(joiner)),
        output_shapes=(
            tf.TensorShape(
                [width, height, feature_builder.get_number_layers()]
            ),  # image array
            tf.TensorShape([]),  # class label
        ),
        output_types=(tf.float64, tf.int64),
    )

    return dataset


def get_generator(stream, length):
    """Return a "batch" of joined (positive and negative) examples.
    The "batch" size should be set to the length of the dataset, in order
    to return a generator of the entire dataset, which is then fed into
    a tf.data.DataSet.
    """

    def generator():
        for el in stream.get_batch(length):
            yield el

    return generator


# def PaddedBatchGenerator3(posStream, negStream, featureBuilder, negRatio, batchSize, pep1Range, pep2Range):
#     """ For positive and negative dataset """
#     width = pep1Range[1]
#     height = pep2Range[1]
#
#     def transform(stream, label):
#         sizeFilter = SizeFilter(stream, pep1Range, pep2Range, hasLabel=False)
#         labeled = Labeler(sizeFilter, label)
#
#         sampler = Sampler(labeled)
#         images = ImageGenerator(sampler, featureBuilder)
#         padding = ImagePadding(images, width, height, padValue=0)
#         return padding
#
#     posImages = transform(posStream, 1)
#     negImages = transform(negStream, 0)
#
#     joiner = Joiner(posImages, negImages, 1-negRatio)
#     batchGenerator = BatchGenerator(joiner, batchSize)
#
#     return batchGenerator


# def padded_batch_generator2(
#     pos_stream,
#     neg_stream,
#     feature_builder,
#     neg_ratio,
#     batch_size,
#     pep1_range,
#     pep2_range,
#     swap=False,
# ):
#     width = pep1_range[1]
#     height = pep2_range[1]

#     pos_size_filter = SizeFilter(pos_stream, pep1_range, pep2_range, has_label=True)
#     pos_filtered1, pos_filtered2 = tee(pos_size_filter)
#     pos_sampler = Sampler(pos_filtered1)

#     if swap:
#         pos_sampler = Swapper(pos_sampler)

#     pos_images = ImageGenerator(pos_sampler, feature_builder)
#     pos_padding = ImagePadding(pos_images, width, height, pad_value=0)

#     neg_size_filter = SizeFilter(neg_stream, pep1_range, pep2_range, has_label=True)
#     neg_filter = PositiveFilter(
#         neg_size_filter, positive_items=pos_filtered2, has_label=True, symmetric=True
#     )
#     neg_sampler = Sampler(neg_filter)

#     if swap:
#         neg_sampler = Swapper(neg_sampler)

#     neg_images = ImageGenerator(neg_sampler, feature_builder)
#     neg_padding = ImagePadding(neg_images, width, height, pad_value=0)

#     joiner = Joiner(pos_padding, neg_padding, 1 - neg_ratio)
#     batch_generator = BatchGenerator(joiner, batch_size, multiple_input=False)

#     return batch_generator
