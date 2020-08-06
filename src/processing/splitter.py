from typing import Tuple

from sklearn.model_selection import train_test_split

from src.data.data_source import DataSource
from src.processing.data_stream import DataStream


def splitter(
    data_source: DataSource, test_size: float, shuffle=True
) -> Tuple[DataStream, DataStream]:
    """ Split a DataSource into a train and test set.

    Needs to be a function to return two streams.

    Parameters
    ----------
    data_source : DataSource
        The input data to divide into a train and test set.
    test_size : float
        The proportion of the dataset to include in the test split.
    shuffle : bool, optional
        Whether the data should be shuffled before splitting, by default True

    Returns
    -------
    (DataStream, DataStream)
        A train and test DataStream.
    """
    data = list(data_source)
    part1, part2 = train_test_split(data, test_size=test_size, shuffle=shuffle)

    return DataStream(part1), DataStream(part2)
