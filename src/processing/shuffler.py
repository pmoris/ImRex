import random

from src.processing.stream import Stream


class Shuffler(Stream):
    def __init__(self, data):
        super().__init__()
        self.data = list(data)
        random.shuffle(self.data)
        self.stream = iter(self.data)

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return iter(self.data)

    def get(self, *args, **kwargs):
        return next(self.stream)

    def on_event(self, event):
        print(event)
        if event == "epoch_end":
            random.shuffle(self.data)
            self.stream = iter(self.data)
