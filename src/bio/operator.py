from functools import lru_cache
import numpy as np

from src.bio.util import scale_matrix


class Operator(object):
    def __repr__(self):
        return self.__class__.__name__ + "()"

    @lru_cache()
    def min_op(self, peptideFeature):
        raise NotImplementedError()

    @lru_cache()
    def max_op(self, peptideFeature):
        raise NotImplementedError()

    def matrix(self, v1, v2):
        raise NotImplementedError()

    def scaledMatrix(self, matrix, upper, peptideFeature):
        return scale_matrix(
            matrix,
            self.min_op(peptideFeature),
            self.max_op(peptideFeature),
            upper=upper,
        )

    def image_matrix(self, v1, v2, peptideFeature):
        m = self.matrix(v1, v2)
        return self.scaledMatrix(m, 255.0, peptideFeature)

    def norm_matrix(self, v1, v2, peptideFeature):
        m = self.matrix(v1, v2)
        return self.scaledMatrix(m, 1.0, peptideFeature)

    def getAmountLayers(self):
        return 1


@lru_cache()
class ProductOperator(Operator):
    @lru_cache()
    def min_op(self, peptideFeature):
        """ The minimum if two values are multiplied """
        if (
            peptideFeature.min < 0 and peptideFeature.max > 0
        ):  # if result can be negative, minimum product is lowest * highest
            return peptideFeature.min * peptideFeature.max
        elif (
            peptideFeature.min < 0 and peptideFeature.max < 0
        ):  # if result positive from negatives, min is closest to zero multiplied
            return peptideFeature.max * peptideFeature.max
        else:  # if result positive from positives, minimum product is lowest * lowest
            return peptideFeature.min * peptideFeature.min

    @lru_cache()
    def max_op(self, peptideFeature):
        """ The maximum if two values are multiplied """
        return max(
            peptideFeature.max * peptideFeature.max,
            peptideFeature.min * peptideFeature.min,
        )

    def matrix(self, v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
        """Performs outer product (or matrix multiplication on column and row vector).

        np.outer(v1, v2) is the standard approach, but for small vectors,
        broadcasting should be faster.
        ref: https://stackoverflow.com/questions/46198104/numpy-multiply-arrays-into-matrix-outer-product

        Note that np.newaxis is an alias for None, which would also work.
        Similarly, the second vector is not required to be broadcast into a matrix, but it makes our intent clearer.

        Parameters
        ----------
        v1 : ndarray
            An ndarray of feature values for each amino acid in a sequence.
        v2 : ndarray
            An ndarray of feature values for each amino acid in a sequence.

        Returns
        -------
        ndarray
            The outer product of the two vectors.
            E.g.
            [1,2] x [3,4] = [[3,4], [6,8]]
        """

        m = v1[:, np.newaxis] * v2[np.newaxis, :]
        return m


@lru_cache()
class LayeredOperator(Operator):
    @lru_cache()
    def min_op(self, peptideFeature):
        raise RuntimeError("min_op not defined for LayeredOperator")

    @lru_cache()
    def max_op(self, peptideFeature):
        raise RuntimeError("max_op not defined for LayeredOperator")

    def matrix(self, v1, v2):
        len1, len2 = len(v1), len(v2)
        v1 = v1.reshape(-1, 1)  # make 2d array with column vector
        v2 = v2.reshape(1, -1)  # make 2d array with row vector
        l1 = np.repeat(v1, len2, axis=1)
        l2 = np.repeat(v2, len1, axis=0)
        l3 = ProductOperator().matrix(v1, v2)
        return np.dstack([l1, l2, l3])

    def scaledMatrix(self, matrix, upper, peptideFeature):
        l1 = scale_matrix(
            np.take(matrix, 0, axis=2),
            peptideFeature.min,
            peptideFeature.max,
            upper=upper,
        )
        l2 = scale_matrix(
            np.take(matrix, 1, axis=2),
            peptideFeature.min,
            peptideFeature.max,
            upper=upper,
        )
        l3 = scale_matrix(
            np.take(matrix, 2, axis=2),
            ProductOperator().min_op(peptideFeature),
            ProductOperator().max_op(peptideFeature),
            upper=upper,
        )
        scaled = np.dstack([l1, l2, l3])
        return scaled

    def getAmountLayers(self):
        return 3


@lru_cache()
class AbsDifferenceOperator(Operator):
    @lru_cache()
    def min_op(self, peptideFeature):
        return 0

    @lru_cache()
    def max_op(self, peptideFeature):
        return abs(peptideFeature.max - peptideFeature.min)

    def matrix(self, v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
        """Compute the pairwise absolute difference between every element in two vectors.

        Parameters
        ----------
        v1 : ndarray
            An ndarray of feature values for each amino acid in a sequence.
        v2 : ndarray
            An ndarray of feature values for each amino acid in a sequence.

        Returns
        -------
        ndarray
            A new matrix where every index is the absolute difference
            of the elements in the corresponding indices of the input vectors.
            E.g.
            [1,2] - [3,4] = [[2,3], [1,2]]
        """

        aa = v1[..., np.newaxis] - v2[np.newaxis, ...]
        aa = np.abs(aa)

        return aa


@lru_cache()
class DifferenceOperator(Operator):
    @lru_cache()
    def min_op(self, peptideFeature):
        return peptideFeature.min - peptideFeature.max

    @lru_cache()
    def max_op(self, peptideFeature):
        return peptideFeature.max - peptideFeature.min

    def matrix(self, v1, v2):
        """ Use subtract rather than product """
        aa = v1[..., np.newaxis] - v2[np.newaxis, ...]
        return aa
