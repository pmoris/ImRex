import pandas as pd

from src.config import PROJECT_ROOT
from src.data.control_cdr3_source import ControlCDR3Source
from src.data.vdjdb_source import VdjdbSource
from src.processing.cv_folds import cv_splitter
from src.processing.negative_sampler import sample_epitope_per_cdr3
from src.processing.splitter import splitter


def test_generate_negatives():
    """ Check that the number of negatives equals the number of positive examples. """
    data_source = VdjdbSource(
        filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )
    data_source.add_pos_labels()
    data_source.generate_negatives_via_shuffling(
        full_dataset_path=PROJECT_ROOT / "src/tests/test_vdjdb.csv"
    )
    assert (
        data_source.data.groupby("y").size()[0]
        == data_source.data.groupby("y").size()[1]
    )


def test_sample_pairs_nan():
    data_source = VdjdbSource(
        filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )
    data_source.add_pos_labels()

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

    # generate negatives through shuffling
    data_source.generate_negatives_via_shuffling(
        full_dataset_path=PROJECT_ROOT / "src/tests/test_vdjdb.csv"
    )

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
    data_source.add_pos_labels()

    full_df = (
        pd.read_csv(
            PROJECT_ROOT / "src/tests/test_vdjdb.csv",
            sep=";",
            usecols=["cdr3", "antigen.epitope"],
        )
        .filter(["cdr3", "antigen.epitope"])
        .drop_duplicates()
        .reset_index(drop=True)
    )

    shuffled_pairs = [
        sample_epitope_per_cdr3(
            cdr3=cdr3,
            df=data_source.data,
            full_df=full_df,
            cdr3_column=data_source.headers["cdr3_header"],
            epitope_column=data_source.headers["epitope_header"],
            seed=seed,
        )
        for seed, cdr3 in enumerate(
            data_source.data[data_source.headers["cdr3_header"]]
        )
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
    data_source.add_pos_labels()

    data_source.generate_negatives_via_shuffling(
        full_dataset_path=PROJECT_ROOT / "src/tests/test_vdjdb.csv"
    )

    train, val = splitter(data_source, test_size=0.2)

    print(pd.DataFrame(train)[1].sort_values().unique())
    assert all(pd.DataFrame(train)[1].sort_values().unique() == [0, 1])
    assert all(pd.DataFrame(val)[1].sort_values().unique() == [0, 1])

    iterations = cv_splitter(data_source=data_source, n_folds=2, cv_type="kfold")
    for train, val in iterations:
        assert all(pd.DataFrame(train)[1].sort_values().unique() == [0, 1])
        assert all(pd.DataFrame(val)[1].sort_values().unique() == [0, 1])

    iterations = cv_splitter(
        data_source=data_source, n_folds=2, cv_type="epitope_grouped",
    )
    for train, val in iterations:
        assert all(pd.DataFrame(train)[1].sort_values().unique() == [0, 1])
        assert all(pd.DataFrame(val)[1].sort_values().unique() == [0, 1])


def test_generate_negatives_from_ref():
    """ Check that the number of negatives equals the number of positive examples. """
    data_source = VdjdbSource(
        filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )
    data_source.add_pos_labels()

    negative_source = ControlCDR3Source(
        filepath=PROJECT_ROOT / "src/tests/test_CDR3_control.tsv",
        min_length=10,
        max_length=20,
    )

    data_source.generate_negatives_from_ref(negative_source)

    assert (
        data_source.data.groupby("y").size()[0]
        == data_source.data.groupby("y").size()[1]
    )


def test_sample_pairs_to_do_neg_ref():
    """ Make sure that the to_do_df only contains negative examples.

    This should be ensured by using keep="first" when removing duplicates.
    """
    data_source = VdjdbSource(
        filepath=PROJECT_ROOT / "src/tests/test_vdjdb.csv",
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )
    data_source.add_pos_labels()

    negative_source = ControlCDR3Source(
        filepath=PROJECT_ROOT / "src/tests/test_CDR3_control.tsv",
        min_length=10,
        max_length=20,
    )

    # sample required number of CDR3 sequences
    amount = data_source.data.shape[0]
    negative_cdr3_series = (
        negative_source.data[negative_source.headers["cdr3_header"]]
        .sample(n=amount)
        .reset_index(drop=True)
        .rename(data_source.headers["cdr3_header"])
    )

    # match with positive epitopes
    negative_df = pd.concat(
        [negative_cdr3_series, data_source.data[data_source.headers["epitope_header"]]],
        axis=1,
    )

    # add class labels
    negative_df["y"] = 0

    # merge with positive dataset
    data_source.data = data_source.data.append(negative_df).reset_index(drop=True)

    # remove false negatives
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
