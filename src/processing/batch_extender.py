import random

from .stream import BatchStream


class BatchExtender(BatchStream):
    def __init__(self, baseStream, extendStream, baseRatio):
        super().__init__(baseStream, extendStream)
        assert 0 <= baseRatio <= 1
        self.baseStream = baseStream
        self.extendStream = extendStream
        self.baseRatio = baseRatio

    def __len__(self):
        return int(len(self.baseStream) // self.baseRatio)

    def getBatch(self, batchSize, *args, **kwargs):
        baseAmount = int(batchSize * self.baseRatio)
        extendAmount = int(batchSize - baseAmount)

        base = self.baseStream.getBatch(baseAmount)

        shape = base[0][0].shape[
            :-1
        ]  # first element = [X, y], we select shape of X (without channels)
        ext = self.extendStream.getBatch(extendAmount, shape=shape)

        all = list()
        all.extend(base)
        all.extend(ext)
        random.shuffle(all)  # shuffle in place
        return all
