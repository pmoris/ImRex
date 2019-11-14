from .data_source_generated import GeneratedDataSource

import random

from bio.peptide import Peptide
from bio.peptide_feature import Hydrophobicity


LENGTH = 10


class DataSourceFeatureGenerated(GeneratedDataSource):
    def generateSample(self):
        # r = random.random()
        # isMatch = 1 if r > 0.5 else 0
        isMatch = random.randint(0, 1)
        x1 = Peptide.random(LENGTH)
        if isMatch:
            x2 = Peptide.generateMatch(x1, Hydrophobicity(), length=10, random_inverse=False)
        else:
            x2 = Peptide.random(LENGTH)

        # return Hydrophobicity().image_matrix(x1, x2).reshape((LENGTH, LENGTH, 1)), isMatch
        return Hydrophobicity().tensor(x1, x2), isMatch

