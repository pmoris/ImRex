import random

from src.processing.stream import BatchStream


class ShapedBatchSampler(BatchStream):
    def __init__(self, *streams, checkStream=None):
        super().__init__(*streams)
        self.streams = streams
        self.checkStream = checkStream

    def __len__(self):
        return None

    def getBatch(self, batchSize, *args, shape=None, **kwargs):
        assert len(shape) == len(self.streams)
        groups = list()
        for size, stream in zip(shape, self.streams):
            groups.append(stream.getGroup(size))

        batch = list()
        while len(batch) < batchSize:
            item = [random.choice(group) for group in groups]
            if self.checkStream:
                key = tuple(len(i) for i in item)
                positive = set(i for i, label in self.checkStream.getGroup(key))
                if tuple(item) in positive:
                    # print(item, "is a positive sample, continuing search for negative sample.")
                    continue  # item is accidentally a positive sample
            batch.append(item)
        return batch
