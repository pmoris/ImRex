from .stream import Stream, TransformStream
from .data_stream import DataStream


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


def Unzipper(stream):
    return tuple(DataStream(l) for l in zip(*stream))
