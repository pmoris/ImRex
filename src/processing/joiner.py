import random

from src.processing.stream import BatchStream


class Joiner(BatchStream):
    def __init__(self, stream1, stream2, ratio):
        super().__init__(stream1, stream2)
        self.stream1 = stream1
        self.stream2 = stream2
        self.ratio = ratio

    def __len__(self):
        source_lengths = [self.stream1.__len__(), self.stream2.__len__()]
        if source_lengths[0] is None:
            return int(source_lengths[1] // (1 - self.ratio))
        elif source_lengths[1] is None:
            return int(source_lengths[0] // self.ratio)
        else:
            lengths = [
                int(source_lengths[1] // (1 - self.ratio)),
                int(source_lengths[0] // self.ratio),
            ]
            return min(lengths)

    def get_batch(self, batch_size, *args, **kwargs):
        amount1 = int(self.ratio * batch_size)
        amount2 = batch_size - amount1

        joined = list()
        for stream, amount in zip([self.stream1, self.stream2], [amount1, amount2]):
            for _ in range(amount):
                item = stream.get(*args, **kwargs)
                # print(item)
                joined.append(item)

        random.shuffle(joined)
        return joined
