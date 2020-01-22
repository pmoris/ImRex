from src.processing.stream import Stream


class Filter(Stream):
    def __init__(self, stream, has_label=True):
        super().__init__(stream)
        self.stream = stream
        self.has_label = has_label

    def __len__(self):
        return None  # size not needed

    def is_valid(self, item):
        raise NotImplementedError

    def get(self, *args, **kwargs):
        while True:
            element = self.stream.get(*args, **kwargs)
            if self.has_label:
                item = element[0]
            else:
                item = element
            if self.is_valid(item):
                return element


class SizeFilter(Filter):
    def __init__(self, stream, *ranges, has_label=True):
        super().__init__(stream, has_label=has_label)
        self.ranges = ranges

    def is_valid(self, item):
        if len(self.ranges) == 1 and not isinstance(item, list):
            item = (item,)

        for pep, (min_size, max_size) in zip(item, self.ranges):
            if len(pep) > max_size or len(pep) < min_size:
                return False

        return True


class PositiveFilter(Filter):
    def __init__(self, stream, positive_items, has_label=True, symmetric=False):
        super().__init__(stream, has_label=has_label)
        self.positive_items = set(i for i, label in positive_items)
        self.symmetric = symmetric

    def is_valid(self, item):
        if tuple(item) in self.positive_items:
            # print("Found a positive item")
            # print(item)
            return False
        if self.symmetric and tuple(reversed(item)) in self.positive_items:
            # print("Found a reversed positive item")
            # print(item)
            return False
        return True
