import random

from .stream import Stream, BatchStream


class Sampler(Stream):
    def __init__(self, stream, infinite=False):
        super().__init__(stream)
        self.samples = list(stream)  # materialize
        self.infinite = infinite

    def __len__(self):
        if self.infinite:
            return None
        return len(self.samples)

    def get(self, *args, **kwargs):
        return random.choice(self.samples)


class GroupSampler(Stream):
    def __init__(self, stream):
        super().__init__(stream)
        self.stream = stream
        self.bins = stream.getGroups()  # materialize

    def __len__(self):
        return self.stream.__len__()

    def get(self, *args, **kwargs):
        keys, weights = zip(*[(k, len(l)) for k, l in self.bins.items()])
        binNr = random.choices(keys, weights=weights, k=1)[0]
        return self.bins[binNr]


class BatchSampler(BatchStream):
    def __init__(self, stream):
        super().__init__(stream)
        self.stream = stream

    def __len__(self):
        return self.stream.__len__()

    def getBatch(self, batchSize, *args, **kwargs):
        items = self.stream.get(*args, **kwargs)

        if len(items) <= batchSize:
            return items

        batch = random.sample(items, batchSize)
        return batch
