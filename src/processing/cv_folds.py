import logging
from typing import Generator, Tuple

import numpy as np
from sklearn.model_selection import GroupKFold, KFold, LeaveOneGroupOut

from src.data.data_source import DataSource
from src.processing.data_stream import DataStream


def cv_splitter(
    data_source: DataSource,
    n_folds: int,
    epitope_grouped: bool = False,
    one_out: bool = False,
    shuffle: bool = True,
) -> Generator[Tuple[DataStream, DataStream], None, None]:
    """Take a DataSource object, split it into folds and yield a tuple of train and test sets as DataStream objects.

    Depending on the provided type, the data will be split into k-folds (after shuffling) or
    split it into k-folds with non-overlapping epitope groups.

    Note that because this it returns a generator object, the logs will only be created upon
    iterating over it, rather than when this function is called the first time.

    Parameters
    ----------
    data_source : DataSource
        The data to be split into folds, contains tuples of sequences tuples and class label 0/1.
    n_folds : int
        The number of folds to create.
    epitope_grouped : bool
        Whether to perform a normal k-fold (False) or grouped k-fold (True) cross-validation, by default False.
    one_out : bool
        Whether to perform a leave-one-epitope out cross-validation, by default False. Requires epitope_grouped to be True.
    shuffle : bool, optional
        Whether to shuffle the data before splitting it into folds, by default True.

    Yields
    ------
    Tuple
        A tuple of DataStream objects, for the train set and for the test set.
    """
    logger = logging.getLogger(__name__)

    # make data indexable
    data_source_array = np.array(list(data_source))

    if epitope_grouped:
        if not one_out:
            cv = GroupKFold(n_splits=n_folds)
            logger.info(
                f"Using a Group K-Fold (k = {n_folds}) cross-validation strategy with the following settings:\n{cv}"
            )
            cv_split = cv.split(
                X=data_source_array,
                groups=data_source.data[data_source.headers["epitope_header"]],
            )
        else:
            cv = LeaveOneGroupOut(n_splits=n_folds)
            logger.info(
                f"Using a Leave One Group Out cross-validation strategy with the following settings:\n{cv}"
            )
            cv_split = cv.split(
                X=data_source_array,
                groups=data_source.data[data_source.headers["epitope_header"]],
            )
    else:
        cv = KFold(n_splits=n_folds, shuffle=shuffle)
        logger.info(
            f"Using a K-Fold (k = {n_folds}) cross-validation strategy with the following settings:\n{cv}"
        )
        cv_split = cv.split(X=data_source_array)

    # loop through folds and return both the indices (for logging purposes) and the associated DataStreams
    for iteration, (train_index, test_index) in enumerate(cv_split):
        train, test = (
            data_source_array[train_index],
            data_source_array[test_index]
            # [data_source_list[i] for i in train_index,
            # [data_source_list[i] for i in test_index]
        )

        cv_log_helper(
            train_index=train_index,
            test_index=test_index,
            iteration=iteration,
            data_source=data_source,
        )

        yield DataStream(train), DataStream(test)


def cv_log_helper(
    train_index: np.ndarray,
    test_index: np.ndarray,
    iteration: int,
    data_source: DataSource,
):
    """Log information about the created folds and optionally exports them.

    Parameters
    ----------
    train_index : np.ndarray
        Indices, corresponding to the vdjdb dataframe, for the train set.
    test_index : np.ndarray
        Indices, corresponding to the vdjdb dataframe, for the test set.
    iteration : int
        Current iteration of cross-validation.
    data_source : DataSource
        The underlying vdjdb data.
    """
    logger = logging.getLogger(__name__)

    # log size of folds and number of labels
    logger.info(
        f"Fold {iteration}: train size: {len(train_index)} -- test size: {len(test_index)}"
    )
    # logger.info(f"Number of unique labels: {np.unique(test[:,1], return_counts=True)}")

    # for non-cdr3 reference type data:
    if "antigen.epitope" in data_source.data.columns:
        # log number of unique epitopes
        train_unique_epitopes = (
            data_source.data.iloc[train_index]
            .loc[:, "antigen.epitope"]
            .unique()
            .shape[0]
        )
        test_unique_epitopes = (
            data_source.data.iloc[train_index]
            .loc[:, "antigen.epitope"]
            .unique()
            .shape[0]
        )
        logger.info(
            f"Unique epitopes in train fold {iteration}: {train_unique_epitopes}"
        )
        logger.info(f"Unique epitopes in test fold {iteration}: {test_unique_epitopes}")

        # log top count epitopes
        train_top_epitopes = (
            data_source.data.iloc[train_index].loc[:, "antigen.epitope"].value_counts()
        )
        logger.info(f"Top epitopes in train fold {iteration}: {train_top_epitopes}")
        test_top_epitopes = (
            data_source.data.iloc[test_index].loc[:, "antigen.epitope"].value_counts()
        )
        logger.info(f"Top epitopes in test fold {iteration}: {test_top_epitopes}")
