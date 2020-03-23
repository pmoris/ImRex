import pandas as pd

from src.config import PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.processing.cv_folds import cv_splitter
from src.processing.splitter import splitter


def test_generate_negatives():
    """ Check that the number of negatives equals the number of positive examples. """
    data_source = VdjdbSource(
        filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )
    data_source.generate_negatives()
    assert (
        data_source.data.groupby("y").size()[0]
        == data_source.data.groupby("y").size()[1]
    )


def test_sample_pairs_nan():
    data_source = VdjdbSource(
        filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )

    # create fake cdr3 sequence that pairs with every epitope in the dataset
    universal_binder_df = pd.DataFrame(
        [
            ("fake_cdr3", epitope, 1)
            for epitope in data_source.data[
                data_source.headers["epitope_header"]
            ].unique()
        ],
        columns=[
            data_source.headers["cdr3_header"],
            data_source.headers["epitope_header"],
            "y",
        ],
    )

    # append as additional positive data
    data_source.data = data_source.data.append(universal_binder_df)

    # generate negatives through sampling
    data_source.generate_negatives()

    # assert that all occurrences of "fake_cdr3" are the original positives, and not new negatives
    fakes = data_source.data.loc[
        data_source.data[data_source.headers["cdr3_header"]] == "fake_cdr3"
    ]
    assert (fakes["y"] == 1).all()  # all fakes are positives
    assert (
        fakes["y"].shape[0]
        == data_source.data[data_source.headers["epitope_header"]].unique().shape[0]
    )  # number of fakes is equal to original number of unique epitopes


def test_sample_pairs_to_do_neg():
    """ Make sure that the to_do_df only contains negative examples.

    The sampling approach should ensure that accidental duplicates of positive pairs (i.e. false negatives)
    never occur, by adding all positive epitope matches to each cdr3 to the list of excluded epitopes when
    sampling a new negative epitope for a given cdr3 sequence.
    """
    data_source = VdjdbSource(
        filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )

    shuffled_pairs = [
        data_source._sample_pairs(cdr3)
        for cdr3 in data_source.data[data_source.headers["cdr3_header"]]
    ]

    # convert list of tuples into dataframe
    shuffled_df = pd.DataFrame(
        shuffled_pairs,
        columns=[
            data_source.headers["cdr3_header"],
            data_source.headers["epitope_header"],
        ],
    )

    # add negative class labels
    shuffled_df["y"] = 0

    # merge with original positive data, ensuring positives are kept at the top of the dataframe
    data_source.data = data_source.data.append(shuffled_df).reset_index(drop=True)

    # extract duplicates
    # NOTE: because the sampling approach ensures that accidental duplicates of positive pairs (i.e. false negatives)
    #       never occur, these will all be accidental duplicate samples of negatives. Therefor, keep="first" is redundant,
    #       but it would result in the positive examples being stored in the dataframe (marking the first (=positives) as True).
    #       This is kept here for historical purposes, because before the epitope was supplied alongside the cdr3 in a zip operation
    #       during sample generation, and the associated positive epitope was used for exclusion purposes.
    to_do_df = data_source.data.loc[
        data_source.data.duplicated(
            subset=[
                data_source.headers["cdr3_header"],
                data_source.headers["epitope_header"],
            ],
            keep="first",
        )
    ]

    assert all(to_do_df["y"] == 0)


def test_cv():
    """ Make sure cv splits contain both negatives and positives. """
    data_source = VdjdbSource(
        filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )

    data_source.generate_negatives()

    train, val = splitter(data_source, test_size=0.2)

    print(pd.DataFrame(train)[1].sort_values().unique())
    assert all(pd.DataFrame(train)[1].sort_values().unique() == [0, 1])
    assert all(pd.DataFrame(val)[1].sort_values().unique() == [0, 1])

    iterations = cv_splitter(
        data_source=data_source, n_folds=2, epitope_grouped=False, one_out=False,
    )
    for train, val in iterations:
        assert all(pd.DataFrame(train)[1].sort_values().unique() == [0, 1])
        assert all(pd.DataFrame(val)[1].sort_values().unique() == [0, 1])

    iterations = cv_splitter(
        data_source=data_source, n_folds=2, epitope_grouped=True, one_out=False,
    )
    for train, val in iterations:
        assert all(pd.DataFrame(train)[1].sort_values().unique() == [0, 1])
        assert all(pd.DataFrame(val)[1].sort_values().unique() == [0, 1])
