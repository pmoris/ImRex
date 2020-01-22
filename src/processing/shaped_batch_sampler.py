import random

from src.processing.stream import BatchStream


class ShapedBatchSampler(BatchStream):
    def __init__(self, *streams, check_stream=None):
        super().__init__(*streams)
        self.streams = streams
        self.check_stream = check_stream

    def __len__(self):
        return None

    def get_batch(self, batch_size, *args, shape=None, **kwargs):
        assert len(shape) == len(self.streams)
        groups = list()
        for size, stream in zip(shape, self.streams):
            groups.append(stream.get_group(size))

        batch = list()
        while len(batch) < batch_size:
            item = [random.choice(group) for group in groups]
            if self.check_stream:
                key = tuple(len(i) for i in item)
                positive = set(i for i, label in self.check_stream.get_group(key))
                if tuple(item) in positive:
                    # print(item, "is a positive sample, continuing search for negative sample.")
                    continue  # item is accidentally a positive sample
            batch.append(item)
        return batch
