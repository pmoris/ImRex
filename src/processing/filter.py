from src.processing.stream import Stream


class Filter(Stream):
    def __init__(self, stream, hasLabel=True):
        super().__init__(stream)
        self.stream = stream
        self.hasLabel = hasLabel

    def __len__(self):
        return None  # size not needed

    def isValid(self, item):
        raise NotImplementedError

    def get(self, *args, **kwargs):
        while True:
            element = self.stream.get(*args, **kwargs)
            if self.hasLabel:
                item = element[0]
            else:
                item = element
            if self.isValid(item):
                return element


class SizeFilter(Filter):
    def __init__(self, stream, *ranges, hasLabel=True):
        super().__init__(stream, hasLabel=hasLabel)
        self.ranges = ranges

    def isValid(self, item):
        if len(self.ranges) == 1 and not isinstance(item, list):
            item = (item,)

        for pep, (minSize, maxSize) in zip(item, self.ranges):
            if len(pep) > maxSize or len(pep) < minSize:
                return False

        return True


class PositiveFilter(Filter):
    def __init__(self, stream, positiveItems, hasLabel=True, symmetric=False):
        super().__init__(stream, hasLabel=hasLabel)
        self.positiveItems = set(i for i, label in positiveItems)
        self.symmetric = symmetric

    def isValid(self, item):
        if tuple(item) in self.positiveItems:
            # print("Found a positive item")
            # print(item)
            return False
        if self.symmetric and tuple(reversed(item)) in self.positiveItems:
            # print("Found a reversed positive item")
            # print(item)
            return False
        return True
