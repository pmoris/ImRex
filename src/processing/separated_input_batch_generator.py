import random

from Bio.SubsMat import MatrixInfo
import numpy as np

from src.definitions.amino_acid_properties import AMINO_ACIDS
from src.processing.batch_generator import BatchGenerator
from src.processing.grouper import GroupedAmountFilter, ShapeGrouper, SizeGrouper
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.sampler import BatchSampler, GroupSampler
from src.processing.shaped_batch_sampler import ShapedBatchSampler
from src.processing.stream import BatchStream, TransformStream
from src.processing.tee import tee
from src.processing.zipper import unzipper

# from src.processing.filter import SizeFilter


def separated_input_batch_generator(
    data_stream,
    neg_ratio,
    batch_size,
    # pep1Range,
    # pep2Range,
    min_amount,
    # negativeStream=None,
):

    # sizeFilter = SizeFilter(dataStream, pep1Range, pep2Range)
    # input1, input2 = Tee(sizeFilter)
    input1, input2 = tee(data_stream)

    shape_grouper = ShapeGrouper(input1)
    grouped_amount_filter = GroupedAmountFilter(shape_grouper, min_amount)
    pos_image_gen = BlossumImageGenerator(grouped_amount_filter)
    group_sampler = GroupSampler(pos_image_gen)
    batch_sampler = BatchSampler(group_sampler)

    label_trimmer = LabelTrimmer(input2)
    peptides1, peptides2 = unzipper(label_trimmer)

    # # random CDR3
    # if negativeStream:
    #     peptides1 = SizeFilter(negativeStream, pep1Range, hasLabel=False)

    peptides1_grouper = SizeGrouper(peptides1, contains_label=False)
    peptides2_grouper = SizeGrouper(peptides2, contains_label=False)
    shaped_batch_sampler = ShapedBatchSampler(
        peptides1_grouper, peptides2_grouper, check_stream=shape_grouper
    )

    negative_labeler = Labeler(shaped_batch_sampler, 0)
    neg_image_gen = BlossumImageGenerator(negative_labeler)

    negative_extender = SeparatedInputBatchExtended(
        batch_sampler, neg_image_gen, 1 - neg_ratio
    )
    batch_generator = BatchGenerator(negative_extender, batch_size, multiple_input=True)

    return batch_generator


class BlossumImageGenerator(TransformStream):
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
        # print(X[0])
        return X, y

    def convert(self, sequence):
        array = np.array([self.features[aa] for aa in sequence])
        # return np.expand_dims(array, -1)
        return array


class SeparatedInputBatchExtended(BatchStream):
    def __init__(self, base_stream, extend_stream, base_ratio):
        super().__init__()
        assert 0 <= base_ratio <= 1
        self.base_stream = base_stream
        self.extend_stream = extend_stream
        self.base_ratio = base_ratio

    def __len__(self):
        return int(len(self.base_stream) // self.base_ratio)

    def get_batch(self, batch_size, *args, **kwargs):
        base_amount = int(batch_size * self.base_ratio)
        extend_amount = int(batch_size - base_amount)

        base = self.base_stream.get_batch(base_amount)

        samples = base[0][0]
        shape = tuple(len(sample) for sample in samples)
        # shape = base[0][0].shape[:-1]         # first element = [X, y], we select shape of X (without channels)
        ext = self.extend_stream.get_batch(extend_amount, shape=shape)

        all_data = list()
        all_data.extend(base)
        all_data.extend(ext)
        random.shuffle(all_data)  # shuffle in place
        return all_data
