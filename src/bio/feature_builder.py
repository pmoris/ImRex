import numpy as np

from src.bio.peptide_feature import ProductOperator


class FeatureBuilder(object):
    def generate_feature(self, source):
        raise NotImplementedError

    def get_number_layers(self):
        raise NotImplementedError


class PeptideFeatureBuilder(FeatureBuilder):
    def __init__(self, operator=ProductOperator()):
        self.operator = operator

    def generate_feature(self, source):
        pep1, pep2 = source
        feature = self.generate_peptides_feature(pep1, pep2)
        return feature

    def generate_peptides_feature(self, pep1, pep2):
        raise NotImplementedError


class SinglePeptideFeatureBuilder(PeptideFeatureBuilder):
    def __init__(self, feature, operator=ProductOperator()):
        super().__init__(operator)
        self.feature = feature

    def generate_peptides_feature(self, pep1, pep2):
        matrix = self.feature.norm_matrix(pep1, pep2, self.operator)
        return np.expand_dims(matrix, axis=-1)

    def get_number_layers(self):
        return self.operator.get_amount_layers()


class CombinedPeptideFeatureBuilder(PeptideFeatureBuilder):
    """Object that contains a list of features and an operator as its attributes.

    Can compute the pairwise feature matrix for two amino acid sequences.
    """

    def __init__(self, features, operator=ProductOperator()):
        super().__init__(operator)
        self.features = features

    def apply_feature(self, feature, pep1, pep2):
        return feature.norm_matrix(pep1, pep2, self.operator)

    def generate_peptides_feature(self, pep1, pep2):
        matrices = []
        for feature in self.features:
            matrix = self.apply_feature(feature, pep1, pep2)
            matrices.append(matrix)
        return np.dstack(matrices)

    def get_number_layers(self):
        if self.operator == "best":
            return sum(f.get_best_operator().get_amount_layers() for f in self.features)
        else:
            return len(self.features) * self.operator.get_amount_layers()
