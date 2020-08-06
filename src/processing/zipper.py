from src.processing.data_stream import DataStream
from src.processing.stream import Stream  # , TransformStream


class Zipper(Stream):
    def __init__(self, *streams):
        super().__init__(*streams)
        self.streams = streams

    def __len__(self):
        lengths = [s.__len__() for s in self.streams]
        if None in lengths:
            return None
        return min(lengths)

    def get(self, *args, **kwargs):
        items = list()
        for stream in self.streams:
            items.append(stream.get(*args, **kwargs))

        return items


def unzipper(stream: Stream):
    """ Split Stream object into separate DataStreams for each element making up a single obervation.

    Input can be a DataStream, TransformStream, etc.

    E.g.

         [
            ('CASGSGAEAFF', 'GILGFVFTL'),
            ('CASSPRDRPLEQYF', 'ELAGIGILTV'),
            ...
         ]  # noqa: RST203

    =>

        (
            ['CASGSGAEAFF','CASSPRDRPLEQYF',...],
            ['GILGFVFTL', 'ELAGIGILTV',...]
        )   # noqa: RST203

    Parameters
    ----------
    stream : Stream
        The Stream to unzip.

    Returns
    -------
    Tuple
        A tuple of DataStreams.
    """
    return tuple(DataStream(l) for l in zip(*stream))
