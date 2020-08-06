from src.processing.batch_extender import BatchExtender
from src.processing.batch_generator import BatchGenerator
from src.processing.grouper import GroupedAmountFilter, ShapeGrouper, SizeGrouper
from src.processing.image_generator import ImageGenerator
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.sampler import BatchSampler, GroupSampler
from src.processing.shaped_batch_sampler import ShapedBatchSampler
from src.processing.tee import tee
from src.processing.zipper import unzipper


def grouped_batch_generator(
    data_stream, feature_builder, neg_ratio, batch_size, min_amount
):
    input1, input2 = tee(data_stream)

    shape_grouper = ShapeGrouper(input1)
    grouped_amount_filter = GroupedAmountFilter(shape_grouper, min_amount)
    pos_image_gen = ImageGenerator(grouped_amount_filter, feature_builder)
    group_sampler = GroupSampler(pos_image_gen)
    batch_sampler = BatchSampler(group_sampler)

    label_trimmer = LabelTrimmer(input2)
    peptides1, peptides2 = unzipper(label_trimmer)

    peptides1_grouper = SizeGrouper(peptides1, contains_label=False)
    peptides2_grouper = SizeGrouper(peptides2, contains_label=False)
    shaped_batch_sampler = ShapedBatchSampler(
        peptides1_grouper, peptides2_grouper, check_stream=shape_grouper
    )

    negative_labeler = Labeler(shaped_batch_sampler, 0)
    neg_image_gen = ImageGenerator(negative_labeler, feature_builder)

    negative_extender = BatchExtender(batch_sampler, neg_image_gen, 1 - neg_ratio)
    batch_generator = BatchGenerator(negative_extender, batch_size)

    return batch_generator


def grouped_batch_generator2(
    data_stream, neg_stream, feature_builder, neg_ratio, batch_size, min_amount
):
    input1, input2 = tee(data_stream)

    shape_grouper = ShapeGrouper(input1)
    grouped_amount_filter = GroupedAmountFilter(shape_grouper, min_amount)
    pos_image_gen = ImageGenerator(grouped_amount_filter, feature_builder)
    group_sampler = GroupSampler(pos_image_gen)
    batch_sampler = BatchSampler(group_sampler)

    label_trimmer = LabelTrimmer(input2)
    cdr3, epitope = unzipper(label_trimmer)

    peptides1_grouper = SizeGrouper(neg_stream, contains_label=False)
    peptides2_grouper = SizeGrouper(epitope, contains_label=False)
    shaped_batch_sampler = ShapedBatchSampler(
        peptides1_grouper, peptides2_grouper, check_stream=shape_grouper
    )

    negative_labeler = Labeler(shaped_batch_sampler, 0)
    neg_image_gen = ImageGenerator(negative_labeler, feature_builder)

    negative_extender = BatchExtender(batch_sampler, neg_image_gen, 1 - neg_ratio)
    batch_generator = BatchGenerator(negative_extender, batch_size)

    return batch_generator
