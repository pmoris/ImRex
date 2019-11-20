from src.processing.batch_generator import BatchGenerator
from src.processing.grouper import ShapeGrouper, GroupedAmountFilter, SizeGrouper
from src.processing.image_generator import ImageGenerator
from src.processing.sampler import GroupSampler, BatchSampler
from src.processing.shaped_batch_sampler import ShapedBatchSampler
from src.processing.batch_extender import BatchExtender
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.tee import Tee
from src.processing.zipper import Unzipper


def GroupedBatchGenerator(dataStream, featureBuilder, negRatio, batchSize, minAmount):
    input1, input2 = Tee(dataStream)

    shapeGrouper = ShapeGrouper(input1)
    groupedAmountFilter = GroupedAmountFilter(shapeGrouper, minAmount)
    posImageGen = ImageGenerator(groupedAmountFilter, featureBuilder)
    groupSampler = GroupSampler(posImageGen)
    batchSampler = BatchSampler(groupSampler)

    labelTrimmer = LabelTrimmer(input2)
    peptides1, peptides2 = Unzipper(labelTrimmer)

    peptides1Grouper = SizeGrouper(peptides1, containsLabel=False)
    peptides2Grouper = SizeGrouper(peptides2, containsLabel=False)
    shapedBatchSampler = ShapedBatchSampler(
        peptides1Grouper, peptides2Grouper, checkStream=shapeGrouper
    )

    negativeLabeler = Labeler(shapedBatchSampler, 0)
    negImageGen = ImageGenerator(negativeLabeler, featureBuilder)

    negativeExtender = BatchExtender(batchSampler, negImageGen, 1 - negRatio)
    batchGenerator = BatchGenerator(negativeExtender, batchSize)

    return batchGenerator


def GroupedBatchGenerator2(
    dataStream, negStream, featureBuilder, negRatio, batchSize, minAmount
):
    input1, input2 = Tee(dataStream)

    shapeGrouper = ShapeGrouper(input1)
    groupedAmountFilter = GroupedAmountFilter(shapeGrouper, minAmount)
    posImageGen = ImageGenerator(groupedAmountFilter, featureBuilder)
    groupSampler = GroupSampler(posImageGen)
    batchSampler = BatchSampler(groupSampler)

    labelTrimmer = LabelTrimmer(input2)
    cdr3, epitope = Unzipper(labelTrimmer)

    peptides1Grouper = SizeGrouper(negStream, containsLabel=False)
    peptides2Grouper = SizeGrouper(epitope, containsLabel=False)
    shapedBatchSampler = ShapedBatchSampler(
        peptides1Grouper, peptides2Grouper, checkStream=shapeGrouper
    )

    negativeLabeler = Labeler(shapedBatchSampler, 0)
    negImageGen = ImageGenerator(negativeLabeler, featureBuilder)

    negativeExtender = BatchExtender(batchSampler, negImageGen, 1 - negRatio)
    batchGenerator = BatchGenerator(negativeExtender, batchSize)

    return batchGenerator
