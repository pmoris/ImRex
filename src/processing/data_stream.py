from src.processing.stream import Stream


class DataStream(Stream):
    """Stream of sequence pair elements, and optionally a class label.

    Format: (('CSARDRTGNTIYF', 'GLCTLVAML'), 1)

    Can be cast to a list or iterated over.

    When piped into Filter or Transform Streams, the underlying
    stream object will get exhausted (via get/next calls) and cannot be re-used
    (unless it is assigned to a new DataStream).
    Note that list or manual iteration will still work.
    """

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
