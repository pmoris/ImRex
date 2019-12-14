import random

from Bio.SubsMat import MatrixInfo
import numpy as np

from src.processing.batch_generator import BatchGenerator
from src.definitions.amino_acid_properties import AMINO_ACIDS
from src.processing.grouper import ShapeGrouper, GroupedAmountFilter, SizeGrouper
from src.processing.sampler import GroupSampler, BatchSampler
from src.processing.shaped_batch_sampler import ShapedBatchSampler
from src.processing.labeler import Labeler, LabelTrimmer
from src.processing.tee import Tee
from src.processing.stream import TransformStream, BatchStream
from src.processing.zipper import Unzipper




def NetTCRBatchGenerator(dataStream, negRatio, batchSize, minAmount):
    input1, input2 = Tee(dataStream)

    shapeGrouper = ShapeGrouper(input1)
    groupedAmountFilter = GroupedAmountFilter(shapeGrouper, minAmount)
    posImageGen = BlossumImageGenerator(groupedAmountFilter)
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
    negImageGen = BlossumImageGenerator(negativeLabeler)

    negativeExtender = NetTCRBatchExtender(batchSampler, negImageGen, 1 - negRatio)
    batchGenerator = BatchGenerator(negativeExtender, batchSize, multipleInput=True)

    return batchGenerator


class BlossumImageGenerator(TransformStream):
    def __init__(self, stream):
        super().__init__(stream)
        self.features = dict()
        for aa in AMINO_ACIDS:
            self.features[aa] = [self._getMatrixEntry(aa, x) for x in AMINO_ACIDS]

    def _getMatrixEntry(self, aa1, aa2):
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


class NetTCRBatchExtender(BatchStream):
    def __init__(self, baseStream, extendStream, baseRatio):
        super().__init__()
        assert 0 <= baseRatio <= 1
        self.baseStream = baseStream
        self.extendStream = extendStream
        self.baseRatio = baseRatio

    def __len__(self):
        return int(len(self.baseStream) // self.baseRatio)

    def getBatch(self, batchSize, *args, **kwargs):
        baseAmount = int(batchSize * self.baseRatio)
        extendAmount = int(batchSize - baseAmount)

        base = self.baseStream.getBatch(baseAmount)

        samples = base[0][0]
        shape = tuple(len(sample) for sample in samples)
        # shape = base[0][0].shape[:-1]         # first element = [X, y], we select shape of X (without channels)
        ext = self.extendStream.getBatch(extendAmount, shape=shape)

        all = list()
        all.extend(base)
        all.extend(ext)
        random.shuffle(all)  # shuffle in place
        return all
