from .joiner import Joiner
from .batch_generator import BatchGenerator
from .sampler import Sampler
from .zipper import Zipper, Unzipper
from .swapper import Swapper
from .labeler import Labeler, LabelTrimmer
from .filter import SizeFilter, PositiveFilter
from .tee import Tee
from .stream import TransformStream

import numpy as np
import keras

from Bio.Alphabet import IUPAC

AMINO_ACIDS = IUPAC.IUPACProtein.letters


def PPILitGenerator(
    dataStream, negRatio, batchSize, pep1Range, pep2Range, swap=False, symmetric=True
):
    width = pep1Range[1]
    height = pep2Range[1]

    sizeFilter = SizeFilter(dataStream, pep1Range, pep2Range)

    input1, input2, input3 = Tee(sizeFilter, amount=3)

    posSampler = Sampler(input1)

    if swap:
        posSampler = Swapper(posSampler)

    posTokenizer = Tokenizer(posSampler)
    posSequencePadder = SequencePadding(posTokenizer, width, height)

    labelTrimmer = LabelTrimmer(input2)
    peptides1, peptides2 = Unzipper(labelTrimmer)

    sampler1 = Sampler(peptides1, infinite=True)
    sampler2 = Sampler(peptides2, infinite=True)
    zipper = Zipper(sampler1, sampler2)
    posFilter = PositiveFilter(
        zipper, positiveItems=input3, hasLabel=False, symmetric=symmetric
    )

    negativeLabeler = Labeler(posFilter, 0)

    if swap:
        negativeLabeler = Swapper(negativeLabeler)

    negTokenizer = Tokenizer(negativeLabeler)
    negSequencePadder = SequencePadding(negTokenizer, width, height)

    joiner = Joiner(posSequencePadder, negSequencePadder, 1 - negRatio)
    batchGenerator = BatchGenerator(joiner, batchSize, multipleInput=True)

    return batchGenerator


def PPILitGenerator2(
    posStream, negStream, negRatio, batchSize, pep1Range, pep2Range, swap=False
):
    width = pep1Range[1]
    height = pep2Range[1]

    tokenizerGenerator = TokenizerGenerator()

    posSizeFilter = SizeFilter(posStream, pep1Range, pep2Range, hasLabel=True)
    posFiltered1, posFiltered2 = Tee(posSizeFilter)

    # posFiltered2 does not need to be swapped
    if swap:
        posFiltered1 = Swapper(posFiltered1)

    posTokenizer = tokenizerGenerator.getTokenizer(posFiltered1)
    posSequencePadder = SequencePadding(posTokenizer, width, height)
    posSampler = Sampler(posSequencePadder)

    negSizeFilter = SizeFilter(negStream, pep1Range, pep2Range, hasLabel=True)
    negFilter = PositiveFilter(
        negSizeFilter, positiveItems=posFiltered2, hasLabel=True, symmetric=True
    )

    if swap:
        negFilter = Swapper(negFilter)

    negTokenizer = tokenizerGenerator.getTokenizer(negFilter)
    negSequencePadder = SequencePadding(negTokenizer, width, height)
    negSampler = Sampler(negSequencePadder)

    joiner = Joiner(posSampler, negSampler, 1 - negRatio)
    batchGenerator = BatchGenerator(joiner, batchSize, multipleInput=True)

    return batchGenerator


class TokenizerGenerator(object):
    def __init__(self):
        self.tokenizer = keras.preprocessing.text.Tokenizer(lower=True, char_level=True)
        self.tokenizer.fit_on_texts(AMINO_ACIDS)

    def getTokenizer(self, stream):
        return Tokenizer(stream, state=self.tokenizer)


class Tokenizer(TransformStream):
    def __init__(self, stream, state=None):
        super().__init__(stream)
        if state is None:
            self.tokenizer = keras.preprocessing.text.Tokenizer(
                lower=True, char_level=True
            )
            self.tokenizer.fit_on_texts(AMINO_ACIDS)
        else:
            self.tokenizer = state

    def transform(self, item, *args, **kwargs):
        sequences, label = item
        a = self.tokenizer.texts_to_sequences(sequences)
        return a, label


class SequencePadding(TransformStream):
    def __init__(self, stream, width, height, padValue=0):
        super().__init__(stream)
        self.width = width
        self.height = height
        self.padValue = padValue

    def transform(self, item, *args, **kwargs):
        (seq1, seq2), label = item
        pad1Length = self.width - len(seq1)
        pad2Length = self.height - len(seq2)

        # seq1 = tuple(seq1) + (self.padValue,) * pad1Length
        # seq2 = tuple(seq2) + (self.padValue,) * pad2Length

        seq1 = (self.padValue,) * pad1Length + tuple(seq1)
        seq2 = (self.padValue,) * pad2Length + tuple(seq2)

        return (seq1, seq2), label
