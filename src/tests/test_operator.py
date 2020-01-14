import numpy as np

from src.bio import operator


def test_product_operator():

    vector_1 = np.array([1, 2, 3])
    vector_2 = np.array([4, 5, 6])

    product_operator = operator.ProductOperator()

    outer_product = product_operator.matrix(vector_1, vector_2)

    expected = np.array([[4, 5, 6], [8, 10, 12], [12, 15, 18]])

    np.testing.assert_array_equal(expected, outer_product)


def test_abs_diff_operator():

    vector_1 = np.array([1, 2, 3])
    vector_2 = np.array([4, 5, 6])

    abs_diff_operator = operator.AbsDifferenceOperator()

    abs_diff = abs_diff_operator.matrix(vector_1, vector_2)

    expected = np.array([[3, 4, 5], [2, 3, 4], [1, 2, 3]])

    np.testing.assert_array_equal(expected, abs_diff)
