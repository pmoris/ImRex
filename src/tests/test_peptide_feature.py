import numpy as np

from src.bio import peptide_feature
from src.definitions.amino_acid_properties import AMINO_ACIDS


def test_features_map():
    assert peptide_feature.features_map
    assert len(peptide_feature.features_map) == 18


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
        [peptide_feature.features_map["charge"]._calculate(aa) for aa in AMINO_ACIDS],
        charge_per_aa,
    )

    # values
    np.testing.assert_almost_equal(
        [*peptide_feature.features_map["charge"].values.values()], charge_per_aa
    )
    np.testing.assert_almost_equal(
        [peptide_feature.features_map["charge"].values.get(aa) for aa in AMINO_ACIDS],
        charge_per_aa,
    )

    # value
    # np.testing.assert_almost_equal(
    #     [peptide_feature.featuresMap["charge"].value(aa) for aa in AMINO_ACIDS],
    #     charge_per_aa,
    # )

    # calculate
    np.testing.assert_almost_equal(
        peptide_feature.features_map["charge"].calculate(AMINO_ACIDS), charge_per_aa
    )

    # calculate manually
    # np.testing.assert_almost_equal(
    #     [peptide_feature.featuresMap["charge"].value(aa) for aa in AMINO_ACIDS],
    #     charge_per_aa,
    # )
    np.testing.assert_almost_equal(
        [peptide_feature.features_map["charge"]._calculate(aa) for aa in AMINO_ACIDS],
        charge_per_aa,
    )
    np.testing.assert_almost_equal(
        [
            peptide_feature.features_map["charge"].values.get(aa, 0)
            for aa in AMINO_ACIDS
        ],
        charge_per_aa,
    )


features_list = [
    peptide_feature.Charge(),
    peptide_feature.Hydrophobicity(),
    peptide_feature.Hydrophilicity(),
    peptide_feature.IsoelectricPoint(),
    peptide_feature.Mass(),
    peptide_feature.Surface(),
    peptide_feature.Flexibility(),
    peptide_feature.Transfer(),
    #  peptide_feature.Prime(),
    peptide_feature.TCRexBasicity(),
    peptide_feature.TCRexHelicity(),
    peptide_feature.TCRexHydrophobicity(),
    peptide_feature.TCRexMutationStability(),
    peptide_feature.AtchleyFactor1(),
    peptide_feature.AtchleyFactor2(),
    peptide_feature.AtchleyFactor3(),
    peptide_feature.AtchleyFactor4(),
    peptide_feature.AtchleyFactor5(),
]


def test__calculate():
    assert all(
        [isinstance(feature._calculate("A"), (int, float)) for feature in features_list]
    )


def test_values():
    for feature in features_list:
        assert len(feature.values) == 20 or len(feature.values) == 23
        assert type(feature.values) == dict


def test_calculate():
    for feature in features_list:
        result = feature.calculate("AAAA")
        assert type(result) == np.ndarray
        assert all([type(i) == np.float64 for i in result])


def test_matrix():
    for feature in features_list:
        result = feature.matrix("AAAA", "RRRR")
        assert type(result) == np.ndarray
        assert result.dtype == np.float64


def test_image_matrix():
    for feature in features_list:
        result = feature.image_matrix("AAAA", "RRRR")
        assert type(result) == np.ndarray
        assert result.dtype == np.uint8
        assert np.logical_and(result >= 0, result <= 255).all()


def test_norm_matrix():
    for feature in features_list:
        result = feature.norm_matrix("AAAA", "RRRR")
        assert type(result) == np.ndarray
        assert result.dtype == np.float64
        assert np.logical_and(result >= 0, result <= 1).all()
