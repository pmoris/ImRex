from src.processing.joiner import Joiner
from src.processing.batch_generator import BatchGenerator
from src.processing.image_generator import ImageGenerator
from src.processing.image_padding import ImagePadding
from src.processing.sampler import Sampler
from src.processing.zipper import Zipper, Unzipper
from src.processing.swapper import Swapper
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.filter import SizeFilter, PositiveFilter
from src.processing.tee import Tee
from src.processing.inverse_map import NoOp


def PaddedBatchGenerator(
    dataStream,
    featureBuilder,
    negRatio,
    batchSize,
    pep1Range,
    pep2Range,
    inverseMap=NoOp(),
    negativeStream=None,
    cacheImages=True,
    swap=False,
):
    """ Standard PaddedBatchGenerator """
    width = pep1Range[1]
    height = pep2Range[1]

    sizeFilter = SizeFilter(dataStream, pep1Range, pep2Range)
    input1, input2, input3 = Tee(sizeFilter, amount=3)

    if not cacheImages:
        input1 = Sampler(input1)

    if swap:
        input1 = Swapper(input1)

    input1 = inverseMap.input(input1)
    posImages = ImageGenerator(input1, featureBuilder)
    posOut = ImagePadding(posImages, width, height, padValue=0)
    posOut = inverseMap.output(posOut)

    if cacheImages:
        posOut = Sampler(posOut)

    labelTrimmer = LabelTrimmer(input2)
    cdr3, epitope = Unzipper(labelTrimmer)

    # random CDR3
    if negativeStream:
        cdr3 = SizeFilter(negativeStream, pep1Range, hasLabel=False)

    sampler1 = Sampler(cdr3, infinite=True)
    sampler2 = Sampler(epitope, infinite=True)
    zipper = Zipper(sampler1, sampler2)
    posFilter = PositiveFilter(zipper, positiveItems=input3, hasLabel=False)

    negativeLabeler = Labeler(posFilter, 0)

    if swap:
        negativeLabeler = Swapper(negativeLabeler)

    negativeLabeler = inverseMap.input(negativeLabeler)
    negImageGen = ImageGenerator(negativeLabeler, featureBuilder)
    negPadding = ImagePadding(negImageGen, width, height, padValue=0)
    negPadding = inverseMap.output(negPadding)

    joiner = Joiner(posOut, negPadding, 1 - negRatio)
    batchGenerator = BatchGenerator(joiner, batchSize)

    return batchGenerator


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


def PaddedBatchGenerator2(
    posStream,
    negStream,
    featureBuilder,
    negRatio,
    batchSize,
    pep1Range,
    pep2Range,
    swap=False,
):
    width = pep1Range[1]
    height = pep2Range[1]

    posSizeFilter = SizeFilter(posStream, pep1Range, pep2Range, hasLabel=True)
    posFiltered1, posFiltered2 = Tee(posSizeFilter)
    posSampler = Sampler(posFiltered1)

    if swap:
        posSampler = Swapper(posSampler)

    posImages = ImageGenerator(posSampler, featureBuilder)
    posPadding = ImagePadding(posImages, width, height, padValue=0)

    negSizeFilter = SizeFilter(negStream, pep1Range, pep2Range, hasLabel=True)
    negFilter = PositiveFilter(
        negSizeFilter, positiveItems=posFiltered2, hasLabel=True, symmetric=True
    )
    negSampler = Sampler(negFilter)

    if swap:
        negSampler = Swapper(negSampler)

    negImages = ImageGenerator(negSampler, featureBuilder)
    negPadding = ImagePadding(negImages, width, height, padValue=0)

    joiner = Joiner(posPadding, negPadding, 1 - negRatio)
    batchGenerator = BatchGenerator(joiner, batchSize, multipleInput=False)

    return batchGenerator
