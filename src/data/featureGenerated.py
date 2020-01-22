import random

from src.bio.peptide import Peptide
from src.bio.peptide_feature import Hydrophobicity
from src.data.data_source_generated import GeneratedDataSource


LENGTH = 10


class DataSourceFeatureGenerated(GeneratedDataSource):
    def generate_sample(self):
        # r = random.random()
        # is_match = 1 if r > 0.5 else 0
        is_match = random.randint(0, 1)
        x1 = Peptide.random(LENGTH)
        if is_match:
            x2 = Peptide.generate_match(
                x1, Hydrophobicity(), length=10, random_inverse=False
            )
        else:
            x2 = Peptide.random(LENGTH)

        # return Hydrophobicity().image_matrix(x1, x2).reshape((LENGTH, LENGTH, 1)), is_match
        return Hydrophobicity().tensor(x1, x2), is_match
