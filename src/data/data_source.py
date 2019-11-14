from processing.stream import Stream


class DataSource(Stream):
    def __init__(self):
        super().__init__()
        self.stream = iter(self)

    def __len__(self):
        raise NotImplementedError

    def __iter__(self):
        raise NotImplementedError

    def get(self, *args, **kwargs):
        return next(self.stream)
