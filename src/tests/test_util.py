import numpy as np

from src.bio import operator, peptide_feature, util


def test_scale_matrix():

    m = np.array(
        [
            [3254.85083169, 5884.43653077, 10566.46556787, 10623.51686787],
            [5884.43653077, 10638.45782041, 19103.08619471, 19206.22909471],
            [10566.46556787, 19103.08619471, 34302.70705801, 34487.91695801],
            [10623.51686787, 19206.22909471, 34487.91695801, 34674.12685801],
        ]
    )

    min_value = operator.ProductOperator().min_op(peptide_feature.Mass())
    max_value = operator.ProductOperator().max_op(peptide_feature.Mass())

    scaled_util = util.scale_matrix(m, oldlower=min_value, oldupper=max_value, upper=1)

    m = np.array(
        [
            [3254.85083169, 5884.43653077, 10566.46556787, 10623.51686787],
            [5884.43653077, 10638.45782041, 19103.08619471, 19206.22909471],
            [10566.46556787, 19103.08619471, 34302.70705801, 34487.91695801],
            [10623.51686787, 19206.22909471, 34487.91695801, 34674.12685801],
        ]
    )

    x_std = (m - min_value) / (max_value - min_value)
    x_scaled = x_std * (1 - 0) + 0

    np.testing.assert_array_almost_equal(scaled_util, x_scaled)
