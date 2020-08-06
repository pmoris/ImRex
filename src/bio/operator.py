from functools import lru_cache

import numpy as np

from src.bio.util import scale_matrix


class Operator(object):
    def __repr__(self):
        return self.__class__.__name__ + "()"

    @lru_cache()
    def min_op(self, peptide_feature):
        raise NotImplementedError()

    @lru_cache()
    def max_op(self, peptide_feature):
        raise NotImplementedError()

    def matrix(self, v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
        """Combine two vectors of amino acid properties in a pairwise fashion, based on the operator.

        Implemented by the specific child classes for each operator.

        Parameters
        ----------
        v1 : np.ndarray
            A vector of amino acid properties (obtained via bio.peptide_feature.calculate).
        v2 : np.ndarray
            A vector of amino acid properties (obtained via bio.peptide_feature.calculate).

        Returns
        -------
        ndarray
        """
        raise NotImplementedError()

    def image_matrix(
        self, v1: np.ndarray, v2: np.ndarray, peptide_feature,
    ) -> np.ndarray:
        """Return a scaled pairwise combined amino acid property matrix for the given sequences and amino acid property.

        Elements are scaled between 0 and 255, where the minimum and maximum value
        are defined by the smallest and largest value that the pairwise combination
        of the amino acid property can take among all 20 amino acids
        (i.e. not the minimum and maximum in the matrix under consideration).

        Parameters
        ----------
        v1 : np.ndarray
            A vector of amino acid properties (obtained via bio.peptide_feature.calculate).
        v2 : np.ndarray
            A vector of amino acid properties (obtained via bio.peptide_feature.calculate).
        peptide_feature : bio.peptide_feature.PeptideFeature
            the amino acid property, used to find the min and max possible value.

        Returns
        -------
        ndarray
            The combined and scaled (0-255) matrix with pairwise combinations of the amino acid values.
        """
        matrix = self.matrix(v1, v2)

        return scale_matrix(
            matrix=matrix,
            oldlower=self.min_op(peptide_feature),
            oldupper=self.max_op(peptide_feature),
            upper=255.0,
        )

    def norm_matrix(
        self, v1: np.ndarray, v2: np.ndarray, peptide_feature,
    ) -> np.ndarray:
        """Return a normalized pairwise combined amino acid property matrix for the given sequences and amino acid property.

        Elements are normalized between 0 and 1, where the minimum and maximum value
        are defined by the smallest and largest value that the pairwise combination
        of the amino acid property can take among all 20 amino acids
        (i.e. not the minimum and maximum in the matrix under consideration).

        Parameters
        ----------
        v1 : np.ndarray
            A vector of amino acid properties (obtained via bio.peptide_feature.calculate).
        v2 : np.ndarray
            A vector of amino acid properties (obtained via bio.peptide_feature.calculate).
        peptide_feature : bio.peptide_feature.PeptideFeature
            the amino acid property, used to find the min and max possible value.

        Returns
        -------
        ndarray
            The combined and normalized (0-1) matrix with pairwise combinations of the amino acid values.
        """
        matrix = self.matrix(v1, v2)

        return scale_matrix(
            matrix=matrix,
            oldlower=self.min_op(peptide_feature),
            oldupper=self.max_op(peptide_feature),
            upper=1.0,
        )

    def get_amount_layers(self):
        return 1


@lru_cache()
class ProductOperator(Operator):
    @lru_cache()
    def min_op(self, peptide_feature):
        """ Find the minimum if two values are multiplied. """
        if (
            peptide_feature.min < 0 and peptide_feature.max > 0
        ):  # if result can be negative, minimum product is lowest * highest
            return peptide_feature.min * peptide_feature.max
        elif (
            peptide_feature.min < 0 and peptide_feature.max < 0
        ):  # if result positive from negatives, min is closest to zero multiplied
            return peptide_feature.max * peptide_feature.max
        else:  # if result positive from positives, minimum product is lowest * lowest
            return peptide_feature.min * peptide_feature.min

    @lru_cache()
    def max_op(self, peptide_feature):
        """ Find the maximum if two values are multiplied. """
        return max(
            peptide_feature.max * peptide_feature.max,
            peptide_feature.min * peptide_feature.min,
        )

    def matrix(self, v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
        """ Perform outer product (or matrix multiplication on column and row vector).

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
    def min_op(self, peptide_feature):
        raise RuntimeError("min_op not defined for LayeredOperator")

    @lru_cache()
    def max_op(self, peptide_feature):
        raise RuntimeError("max_op not defined for LayeredOperator")

    def matrix(self, v1, v2):
        len1, len2 = len(v1), len(v2)
        v1 = v1.reshape(-1, 1)  # make 2d array with column vector
        v2 = v2.reshape(1, -1)  # make 2d array with row vector
        l1 = np.repeat(v1, len2, axis=1)
        l2 = np.repeat(v2, len1, axis=0)
        l3 = ProductOperator().matrix(v1, v2)
        return np.dstack([l1, l2, l3])

    def scaled_matrix(self, matrix, upper, peptide_feature):
        l1 = scale_matrix(
            np.take(matrix, 0, axis=2),
            oldlower=peptide_feature.min,
            oldupper=peptide_feature.max,
            upper=upper,
        )
        l2 = scale_matrix(
            np.take(matrix, 1, axis=2),
            oldlower=peptide_feature.min,
            oldupper=peptide_feature.max,
            upper=upper,
        )
        l3 = scale_matrix(
            np.take(matrix, 2, axis=2),
            ProductOperator().min_op(peptide_feature),
            ProductOperator().max_op(peptide_feature),
            upper=upper,
        )
        scaled = np.dstack([l1, l2, l3])
        return scaled

    def get_amount_layers(self):
        return 3


@lru_cache()
class AbsDifferenceOperator(Operator):
    @lru_cache()
    def min_op(self, peptide_feature):
        return 0

    @lru_cache()
    def max_op(self, peptide_feature):
        return abs(peptide_feature.max - peptide_feature.min)

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
    def min_op(self, peptide_feature):
        return peptide_feature.min - peptide_feature.max

    @lru_cache()
    def max_op(self, peptide_feature):
        return peptide_feature.max - peptide_feature.min

    def matrix(self, v1, v2):
        """ Use subtract rather than product. """
        aa = v1[..., np.newaxis] - v2[np.newaxis, ...]
        return aa
