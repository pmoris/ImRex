from src.config import PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.preprocessing.decoy_epitopes import create_decoy_epitopes


def test_decoy_epitopes():
    """Check that all epitopes were converted to a new unique decoy."""

    decoy_df = create_decoy_epitopes(PROJECT_ROOT / "src/tests/test_vdjdb.csv")

    data_source = VdjdbSource(
        filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )

    assert set(decoy_df["antigen.epitope"]).isdisjoint(
        data_source.data["antigen.epitope"]
    )
    assert len(decoy_df["antigen.epitope"].unique()) == len(
        data_source.data["antigen.epitope"].unique()
    )
