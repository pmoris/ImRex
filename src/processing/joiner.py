import random

from .stream import BatchStream


class Joiner(BatchStream):
    def __init__(self, stream1, stream2, ratio):
        super().__init__(stream1, stream2)
        self.stream1 = stream1
        self.stream2 = stream2
        self.ratio = ratio

    def __len__(self):
        sourceLengths = [self.stream1.__len__(), self.stream2.__len__()]
        if sourceLengths[0] is None:
            return int(sourceLengths[1] // (1 - self.ratio))
        elif sourceLengths[1] is None:
            return int(sourceLengths[0] // self.ratio)
        else:
            lengths = [
                int(sourceLengths[1] // (1 - self.ratio)),
                int(sourceLengths[0] // self.ratio),
            ]
            return min(lengths)

    def getBatch(self, batchSize, *args, **kwargs):
        amount1 = int(self.ratio * batchSize)
        amount2 = batchSize - amount1

        joined = list()
        for stream, amount in zip([self.stream1, self.stream2], [amount1, amount2]):
            for _ in range(amount):
                item = stream.get(*args, **kwargs)
                # print(item)
                joined.append(item)

        random.shuffle(joined)
        return joined
