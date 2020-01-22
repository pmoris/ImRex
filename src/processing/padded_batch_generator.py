from src.processing.joiner import Joiner
from src.processing.batch_generator import BatchGenerator
from src.processing.image_generator import ImageGenerator
from src.processing.image_padding import ImagePadding
from src.processing.sampler import Sampler
from src.processing.zipper import Zipper, unzipper
from src.processing.swapper import Swapper
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.filter import SizeFilter, PositiveFilter
from src.processing.tee import tee
from src.processing.inverse_map import NoOp


def padded_batch_generator(
    data_stream,
    feature_builder,
    neg_ratio,
    batch_size,
    pep1_range,
    pep2_range,
    inverse_map=NoOp(),
    negative_stream=None,
    cache_images=True,
    swap=False,
):
    """ Standard PaddedBatchGenerator """
    width = pep1_range[1]
    height = pep2_range[1]

    size_filter = SizeFilter(data_stream, pep1_range, pep2_range)
    input1, input2, input3 = tee(size_filter, amount=3)

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

    label_trimmer = LabelTrimmer(input2)
    cdr3, epitope = unzipper(label_trimmer)

    # random CDR3
    if negative_stream:
        cdr3 = SizeFilter(negative_stream, pep1_range, has_label=False)

    sampler1 = Sampler(cdr3, infinite=True)
    sampler2 = Sampler(epitope, infinite=True)
    zipper = Zipper(sampler1, sampler2)
    pos_filter = PositiveFilter(zipper, positive_items=input3, has_label=False)

    negative_labeler = Labeler(pos_filter, 0)

    if swap:
        negative_labeler = Swapper(negative_labeler)

    negative_labeler = inverse_map.input(negative_labeler)
    neg_image_gen = ImageGenerator(negative_labeler, feature_builder)
    neg_padding = ImagePadding(neg_image_gen, width, height, pad_value=0)
    neg_padding = inverse_map.output(neg_padding)

    joiner = Joiner(pos_out, neg_padding, 1 - neg_ratio)
    batch_generator = BatchGenerator(joiner, batch_size)

    return batch_generator


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


def padded_batch_generator2(
    pos_stream,
    neg_stream,
    feature_builder,
    neg_ratio,
    batch_size,
    pep1_range,
    pep2_range,
    swap=False,
):
    width = pep1_range[1]
    height = pep2_range[1]

    pos_size_filter = SizeFilter(pos_stream, pep1_range, pep2_range, has_label=True)
    pos_filtered1, pos_filtered2 = tee(pos_size_filter)
    pos_sampler = Sampler(pos_filtered1)

    if swap:
        pos_sampler = Swapper(pos_sampler)

    pos_images = ImageGenerator(pos_sampler, feature_builder)
    pos_padding = ImagePadding(pos_images, width, height, pad_value=0)

    neg_size_filter = SizeFilter(neg_stream, pep1_range, pep2_range, has_label=True)
    neg_filter = PositiveFilter(
        neg_size_filter, positive_items=pos_filtered2, has_label=True, symmetric=True
    )
    neg_sampler = Sampler(neg_filter)

    if swap:
        neg_sampler = Swapper(neg_sampler)

    neg_images = ImageGenerator(neg_sampler, feature_builder)
    neg_padding = ImagePadding(neg_images, width, height, pad_value=0)

    joiner = Joiner(pos_padding, neg_padding, 1 - neg_ratio)
    batch_generator = BatchGenerator(joiner, batch_size, multiple_input=False)

    return batch_generator
