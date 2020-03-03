from functools import lru_cache


class _Stream(object):
    def __init__(self, *sources):
        self.sources = sources

    def __len__(self):
        raise NotImplementedError

    def send_event(self, event):
        self.on_event(event)
        for source in self.sources:
            source.send_event(event)

    def on_event(self, event):
        pass


class Stream(_Stream):
    def __init__(self, *sources):
        super().__init__(*sources)

    def __iter__(self):
        return self

    def __next__(self):
        return self.get()

    def get(self, *args, **kwargs):
        raise NotImplementedError


class BatchStream(_Stream):
    def __init__(self, *sources):
        super().__init__(*sources)

    def get_batch(self, batch_size, *args, **kwargs):
        raise NotImplementedError


class GroupedStream(_Stream):
    def __init__(self, *sources):
        super().__init__(*sources)

    def get_group(self, key):
        return self.get_groups().get(key)

    @lru_cache()
    def get_groups(self):
        raise NotImplementedError

    def __iter__(self):
        return ((k, v) for k, v in self.get_groups().items())


class TransformStream(Stream, GroupedStream, BatchStream):
    """ Iterable that transforms each element upon iteration. """

    def __init__(self, stream):
        super().__init__(stream)
        self.stream = stream

    def __len__(self):
        return self.stream.__len__()

    def transform(self, item, *args, **kwargs):
        raise NotImplementedError

    def get(self, *args, **kwargs):
        item = self.stream.get(*args, **kwargs)
        return self.transform(item, *args, **kwargs)

    def get_groups(self):
        # item = [X, y] = [(pep1, pep2), y]
        return {
            k: [(self.transform(item)) for item in bin]
            for k, bin in self.stream.get_groups().items()
        }

    def get_batch(self, batch_size, *args, **kwargs):
        batch = self.stream.get_batch(batch_size, *args, **kwargs)
        return ((self.transform(item)) for item in batch)
