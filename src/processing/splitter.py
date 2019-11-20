from sklearn.model_selection import train_test_split

from src.processing.data_stream import DataStream


def Splitter(dataSource, ratio, shuffle=True):
    """ Needs to be a function to return two streams """
    data = list(dataSource)
    part1, part2 = train_test_split(data, test_size=ratio, shuffle=shuffle)

    return DataStream(part1), DataStream(part2)
