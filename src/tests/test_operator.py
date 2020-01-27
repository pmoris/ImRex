import numpy as np

from src.bio import peptide_feature, operator
from src.definitions.amino_acid_properties import AMINO_ACIDS


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


def test_min_op():
    # product
    np.testing.assert_almost_equal(
        operator.ProductOperator().min_op(peptide_feature.Charge()),
        -0.9995470362780428,
    )
    np.testing.assert_almost_equal(
        operator.ProductOperator().max_op(peptide_feature.Charge()), 1.00314089534477
    )
    np.testing.assert_almost_equal(
        operator.ProductOperator().min_op(peptide_feature.Mass()), 3254.8508316899997
    )
    np.testing.assert_almost_equal(
        operator.ProductOperator().max_op(peptide_feature.Mass()), 34674.12685801
    )

    # abs diff
    np.testing.assert_almost_equal(
        operator.AbsDifferenceOperator().min_op(peptide_feature.Charge()), 0,
    )
    np.testing.assert_almost_equal(
        operator.AbsDifferenceOperator().max_op(peptide_feature.Charge()),
        1.9995502045447453,
    )
    np.testing.assert_almost_equal(
        operator.AbsDifferenceOperator().min_op(peptide_feature.Mass()), 0
    )
    np.testing.assert_almost_equal(
        operator.AbsDifferenceOperator().max_op(peptide_feature.Mass()), 129.1586,
    )


def test_image_matrix():

    # abs diff
    v1 = np.array([57.0513, 103.1429, 186.2099])
    v2 = v1
    m = operator.AbsDifferenceOperator().image_matrix(
        v1=v1, v2=v2, peptide_feature=peptide_feature.Mass()
    )

    np.testing.assert_almost_equal(m.max(), 255.0)
    np.testing.assert_almost_equal(m.min(), 0.0)
    np.testing.assert_array_almost_equal(
        m,
        np.array(
            [
                [0.0, 90.99942242, 255.0],
                [90.99942242, 0.0, 164.00057758],
                [255.0, 164.00057758, 0.0],
            ]
        ),
    )

    # product
    m = operator.ProductOperator().image_matrix(
        v1=v1, v2=v2, peptide_feature=peptide_feature.Mass()
    )

    np.testing.assert_almost_equal(m.max(), 255.0)
    np.testing.assert_almost_equal(m.min(), 0.0)
    np.testing.assert_array_almost_equal(
        m,
        np.array(
            [
                [0.0, 21.34181426, 59.80436461],
                [21.34181426, 59.92562593, 129.46197276],
                [59.80436461, 129.46197276, 255.0],
            ]
        ),
    )

    # all features min/max
    for feature in [
        i for i in peptide_feature.features_map.values() if i.name != "Prime"
    ]:

        v1 = feature.calculate(AMINO_ACIDS)
        v2 = v1
        m = operator.ProductOperator().image_matrix(
            v1=v1, v2=v2, peptide_feature=feature
        )

        print("\n\nMATRIX M\n\n")
        print(m)
        print(v1, v2)

        np.testing.assert_almost_equal(m.max(), 255.0)
        np.testing.assert_almost_equal(m.min(), 0.0)

    # manual calculation
    v1 = np.array([57.0513, 103.1429, 185.2099, 186.2099])
    v2 = v1
    m = operator.ProductOperator().image_matrix(
        v1=v1, v2=v2, peptide_feature=peptide_feature.Mass()
    )

    m_raw = operator.ProductOperator().matrix(v1=v1, v2=v2)

    x_std = (m_raw - operator.ProductOperator().min_op(peptide_feature.Mass())) / (
        operator.ProductOperator().max_op(peptide_feature.Mass())
        - operator.ProductOperator().min_op(peptide_feature.Mass())
    )
    x_scaled = x_std * (255 - 0) + 0

    np.testing.assert_array_almost_equal(m, x_scaled)


def test_norm_matrix():

    # abs diff
    v1 = np.array([57.0513, 103.1429, 186.2099])
    v2 = v1
    m = operator.AbsDifferenceOperator().norm_matrix(
        v1=v1, v2=v2, peptide_feature=peptide_feature.Mass()
    )

    np.testing.assert_almost_equal(m.max(), 1.0)
    np.testing.assert_almost_equal(m.min(), 0.0)
    np.testing.assert_array_almost_equal(
        m,
        np.array(
            [
                [0.0, 0.35686048, 1.0],
                [0.35686048, 0.0, 0.64313952],
                [1.0, 0.64313952, 0.0],
            ]
        ),
    )

    # product
    m = operator.ProductOperator().norm_matrix(
        v1=v1, v2=v2, peptide_feature=peptide_feature.Mass()
    )

    np.testing.assert_almost_equal(m.max(), 1.0)
    np.testing.assert_almost_equal(m.min(), 0.0)
    np.testing.assert_array_almost_equal(
        m,
        np.array(
            [
                [0.0, 0.08369339, 0.23452692],
                [0.08369339, 0.23500245, 0.50769401],
                [0.23452692, 0.50769401, 1.0],
            ]
        ),
    )

    # all features min/max
    for feature in [
        i for i in peptide_feature.features_map.values() if i.name != "Prime"
    ]:
        v1 = feature.calculate(AMINO_ACIDS)
        v2 = v1
        m = operator.ProductOperator().norm_matrix(
            v1=v1, v2=v2, peptide_feature=feature
        )
        np.testing.assert_almost_equal(m.max(), 1.0)
        np.testing.assert_almost_equal(m.min(), 0.0)

    # manual calculation
    v1 = np.array([57.0513, 103.1429, 185.2099, 186.2099])
    v2 = v1

    m = operator.ProductOperator().norm_matrix(
        v1=v1, v2=v2, peptide_feature=peptide_feature.Mass()
    )

    m_raw = operator.ProductOperator().matrix(v1=v1, v2=v2)

    x_std = (m_raw - operator.ProductOperator().min_op(peptide_feature.Mass())) / (
        operator.ProductOperator().max_op(peptide_feature.Mass())
        - operator.ProductOperator().min_op(peptide_feature.Mass())
    )
    x_scaled = x_std * (1 - 0) + 0

    np.testing.assert_array_almost_equal(m, x_scaled)


def test_matrix_type():
    """
    In-place operations do not change the dtype of the container array. Since the desired normalized values are floats, the matrix need to have floating-point point dtype before the in-place operations are performed.

    All peptide_feature.calculate methods, which proceed any Operator().matrix() methods, should return floats by default,
    since all amino acid properties are floats.

    This test simply loops through the different features and checks whether the generated matrices have the dtype

    """
    for feature in [
        i for i in peptide_feature.features_map.values() if i.name != "Prime"
    ]:
        for operator_looper in []:

            v1 = feature.calculate(AMINO_ACIDS)
            v2 = v1
            m = operator.operator_looper.matrix(v1=v1, v2=v2, peptide_feature=feature)
            assert isinstance(m, np.float64)
