from bio.peptide_feature import *

import numpy as np


class FeatureBuilder(object):
    def generateFeature(self, source):
        raise NotImplementedError

    def getNumberLayers(self):
        raise NotImplementedError


class PeptideFeatureBuilder(FeatureBuilder):
    def __init__(self, operator=ProductOperator()):
        self.operator = operator

    def generateFeature(self, source):
        pep1, pep2 = source
        feature = self.generatePeptidesFeature(pep1, pep2)
        return feature

    def generatePeptidesFeature(self, pep1, pep2):
        raise NotImplementedError


class SinglePeptideFeatureBuilder(PeptideFeatureBuilder):
    def __init__(self, feature, operator=ProductOperator()):
        super().__init__(operator)
        self.feature = feature

    def generatePeptidesFeature(self, pep1, pep2):
        # print(pep1, pep2)
        matrix = self.feature.norm_matrix(pep1, pep2, self.operator)
        # print(matrix)
        return np.expand_dims(matrix, axis=-1)

    def getNumberLayers(self):
        return self.operator.getAmountLayers()


class CombinedPeptideFeatureBuilder(PeptideFeatureBuilder):
    def __init__(self, features, operator=ProductOperator()):
        super().__init__(operator)
        self.features = features

    def applyFeature(self, feature, pep1, pep2):
        return feature.norm_matrix(pep1, pep2, self.operator)

    def generatePeptidesFeature(self, pep1, pep2):
        matrices = []
        for feature in self.features:
            matrix = self.applyFeature(feature, pep1, pep2)
            matrices.append(matrix)
        return np.dstack(matrices)

    def getNumberLayers(self):
        if self.operator == 'best':
            return sum(f.getBestOperator().getAmountLayers() for f in self.features)
        else:
            return len(self.features) * self.operator.getAmountLayers()

