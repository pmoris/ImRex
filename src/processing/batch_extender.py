import random

from src.processing.stream import BatchStream


class BatchExtender(BatchStream):
    def __init__(self, base_stream, extend_stream, base_ratio):
        super().__init__(base_stream, extend_stream)
        assert 0 <= base_ratio <= 1
        self.base_stream = base_stream
        self.extend_stream = extend_stream
        self.base_ratio = base_ratio

    def __len__(self):
        return int(len(self.base_stream) // self.base_ratio)

    def get_batch(self, batch_size, *args, **kwargs):
        base_amount = int(batch_size * self.base_ratio)
        extend_amount = int(batch_size - base_amount)

        base = self.base_stream.get_batch(base_amount)

        shape = base[0][0].shape[
            :-1
        ]  # first element = [X, y], we select shape of X (without channels)
        ext = self.extend_stream.get_batch(extend_amount, shape=shape)

        all_data = list()
        all_data.extend(base)
        all_data.extend(ext)
        random.shuffle(all_data)  # shuffle in place
        return all_data
