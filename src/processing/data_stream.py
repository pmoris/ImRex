from src.processing.stream import Stream


class DataStream(Stream):
    def __init__(self, data):
        super().__init__()
        self.data = list(data)
        self.stream = iter(self.data)

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return iter(self.data)

    def get(self, *args, **kwargs):
        return next(self.stream)
