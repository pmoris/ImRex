from functools import lru_cache


class _Stream(object):
    def __init__(self, *sources):
        self.sources = sources

    def __len__(self):
        raise NotImplementedError

    def sendEvent(self, event):
        self.onEvent(event)
        for source in self.sources:
            source.sendEvent(event)

    def onEvent(self, event):
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

    def getBatch(self, batchSize, *args, **kwargs):
        raise NotImplementedError


class GroupedStream(_Stream):
    def __init__(self, *sources):
        super().__init__(*sources)

    def getGroup(self, key):
        return self.getGroups().get(key)

    @lru_cache()
    def getGroups(self):
        raise NotImplementedError

    def __iter__(self):
        return ((k, v) for k, v in self.getGroups().items())


class TransformStream(Stream, GroupedStream, BatchStream):
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

    def getGroups(self):
        # item = [X, y] = [(pep1, pep2), y]
        return {k: [(self.transform(item)) for item in bin] for k, bin in self.stream.getGroups().items()}

    def getBatch(self, batchSize, *args, **kwargs):
        batch = self.stream.getBatch(batchSize, *args, **kwargs)
        return ((self.transform(item)) for item in batch)
