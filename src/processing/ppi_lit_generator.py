from tensorflow import keras

from src.definitions.amino_acid_properties import AMINO_ACIDS
from src.processing.batch_generator import BatchGenerator
from src.processing.filter import PositiveFilter, SizeFilter
from src.processing.joiner import Joiner
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.sampler import Sampler
from src.processing.stream import TransformStream
from src.processing.swapper import Swapper
from src.processing.tee import tee
from src.processing.zipper import unzipper, Zipper


def ppi_lit_generator(
    data_stream,
    neg_ratio,
    batch_size,
    pep1_range,
    pep2_range,
    swap=False,
    symmetric=True,
):
    width = pep1_range[1]
    height = pep2_range[1]

    size_filter = SizeFilter(data_stream, pep1_range, pep2_range)

    input1, input2, input3 = tee(size_filter, amount=3)

    pos_sampler = Sampler(input1)

    if swap:
        pos_sampler = Swapper(pos_sampler)

    pos_tokenizer = Tokenizer(pos_sampler)
    pos_sequence_padder = SequencePadding(pos_tokenizer, width, height)

    label_trimmer = LabelTrimmer(input2)
    peptides1, peptides2 = unzipper(label_trimmer)

    sampler1 = Sampler(peptides1, infinite=True)
    sampler2 = Sampler(peptides2, infinite=True)
    zipper = Zipper(sampler1, sampler2)
    pos_filter = PositiveFilter(
        zipper, positive_items=input3, has_label=False, symmetric=symmetric
    )

    negative_labeler = Labeler(pos_filter, 0)

    if swap:
        negative_labeler = Swapper(negative_labeler)

    neg_tokenizer = Tokenizer(negative_labeler)
    neg_sequence_padder = SequencePadding(neg_tokenizer, width, height)

    joiner = Joiner(pos_sequence_padder, neg_sequence_padder, 1 - neg_ratio)
    batch_generator = BatchGenerator(joiner, batch_size, multiple_input=True)

    return batch_generator


def ppi_lit_generator2(
    pos_stream, neg_stream, neg_ratio, batch_size, pep1_range, pep2_range, swap=False
):
    width = pep1_range[1]
    height = pep2_range[1]

    tokenizer_generator = TokenizerGenerator()

    pos_size_filter = SizeFilter(pos_stream, pep1_range, pep2_range, has_label=True)
    pos_filtered1, pos_filtered2 = tee(pos_size_filter)

    # pos_filtered2 does not need to be swapped
    if swap:
        pos_filtered1 = Swapper(pos_filtered1)

    pos_tokenizer = tokenizer_generator.get_tokenizer(pos_filtered1)
    pos_sequence_padder = SequencePadding(pos_tokenizer, width, height)
    pos_sampler = Sampler(pos_sequence_padder)

    neg_size_filter = SizeFilter(neg_stream, pep1_range, pep2_range, has_label=True)
    neg_filter = PositiveFilter(
        neg_size_filter, positive_items=pos_filtered2, has_label=True, symmetric=True
    )

    if swap:
        neg_filter = Swapper(neg_filter)

    neg_tokenizer = tokenizer_generator.get_tokenizer(neg_filter)
    neg_sequence_padder = SequencePadding(neg_tokenizer, width, height)
    neg_sampler = Sampler(neg_sequence_padder)

    joiner = Joiner(pos_sampler, neg_sampler, 1 - neg_ratio)
    batch_generator = BatchGenerator(joiner, batch_size, multiple_input=True)

    return batch_generator


class TokenizerGenerator(object):
    def __init__(self):
        self.tokenizer = keras.preprocessing.text.Tokenizer(lower=True, char_level=True)
        self.tokenizer.fit_on_texts(AMINO_ACIDS)

    def get_tokenizer(self, stream):
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
    def __init__(self, stream, width, height, pad_value=0):
        super().__init__(stream)
        self.width = width
        self.height = height
        self.pad_value = pad_value

    def transform(self, item, *args, **kwargs):
        (seq1, seq2), label = item
        pad1_length = self.width - len(seq1)
        pad2_length = self.height - len(seq2)

        # seq1 = tuple(seq1) + (self.padValue,) * pad1_length
        # seq2 = tuple(seq2) + (self.padValue,) * pad2_length

        seq1 = (self.pad_value,) * pad1_length + tuple(seq1)
        seq2 = (self.pad_value,) * pad2_length + tuple(seq2)

        return (seq1, seq2), label
