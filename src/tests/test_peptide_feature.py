import numpy as np

from src.bio import peptide_feature
from src.definitions.amino_acid_properties import AMINO_ACIDS


def test_features_map():
    assert peptide_feature.featuresMap
    assert len(peptide_feature.featuresMap) == 18


def test_charge():

    charge_per_aa = [
        -0.0020157006072527572,
        -0.06399041653261084,
        -1.001569216452248,
        -1.000240577861443,
        -0.0020157006072527572,
        -0.0020157006072527572,
        0.08889339030183807,
        -0.0020157006072527572,
        0.9976892655407432,
        -0.0020157006072527572,
        -0.0020157006072527572,
        -0.0020157006072527572,
        -0.0020157006072527572,
        -0.0020157006072527572,
        0.9979809880924972,
        -0.0020157006072527572,
        -0.0020157006072527572,
        -0.0020157006072527572,
        -0.0020157006072527572,
        -0.0028661148255656466,
    ]

    # _calculate
    np.testing.assert_almost_equal(
        [peptide_feature.featuresMap["charge"]._calculate(aa) for aa in AMINO_ACIDS],
        charge_per_aa,
    )

    # values
    np.testing.assert_almost_equal(
        [*peptide_feature.featuresMap["charge"].values.values()], charge_per_aa
    )
    np.testing.assert_almost_equal(
        [peptide_feature.featuresMap["charge"].values.get(aa) for aa in AMINO_ACIDS],
        charge_per_aa,
    )

    # value
    # np.testing.assert_almost_equal(
    #     [peptide_feature.featuresMap["charge"].value(aa) for aa in AMINO_ACIDS],
    #     charge_per_aa,
    # )

    # calculate
    np.testing.assert_almost_equal(
        peptide_feature.featuresMap["charge"].calculate(AMINO_ACIDS), charge_per_aa
    )

    # calculate manually
    # np.testing.assert_almost_equal(
    #     [peptide_feature.featuresMap["charge"].value(aa) for aa in AMINO_ACIDS],
    #     charge_per_aa,
    # )
    np.testing.assert_almost_equal(
        [peptide_feature.featuresMap["charge"]._calculate(aa) for aa in AMINO_ACIDS],
        charge_per_aa,
    )
    np.testing.assert_almost_equal(
        [peptide_feature.featuresMap["charge"].values.get(aa, 0) for aa in AMINO_ACIDS],
        charge_per_aa,
    )
