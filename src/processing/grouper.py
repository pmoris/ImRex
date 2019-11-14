from collections import defaultdict
from functools import lru_cache

from .stream import GroupedStream


class Grouper(GroupedStream):
    def __init__(self, stream):
        super().__init__(stream)
        self.stream = stream

    def __len__(self):
        return len(self.stream)

    def getKey(self, item):
        raise NotImplementedError

    @lru_cache()
    def getGroups(self):
        bins = defaultdict(list)
        for item in self.stream:
            key = self.getKey(item)
            bins[key].append(item)
        return bins


class ShapeGrouper(Grouper):
    def __init__(self, stream, containsLabel=True):
        super().__init__(stream)
        self.containsLabel = containsLabel

    def getKey(self, item):
        if self.containsLabel:
            item = item[0]
        return tuple(len(i) for i in item)


class SizeGrouper(Grouper):
    def __init__(self, stream, containsLabel=True):
        super().__init__(stream)
        self.containsLabel = containsLabel

    def getKey(self, item):
        if self.containsLabel:
            item = item[0]
        return len(item)


class GroupedAmountFilter(GroupedStream):
    def __init__(self, stream, minAmount):
        super().__init__(stream)
        self.stream = stream
        self.minAmount = minAmount

    def __len__(self):
        return sum([len(bin) for bin in self.getGroups().values()])

    @lru_cache()
    def getGroups(self):
        # discarded = 0
        # discardedBins = 0
        # for k, bin in self.stream:
        #     # print(k, "\t", len(bin))
        #     if len(bin) < self.minAmount:
        #         # print("Discard")
        #         discardedBins += 1
        #         discarded += len(bin)
        #
        # print("\nTotal discarded:", discarded)
        # print("Discarded nins:", discardedBins)
        # print()

        return {k: bin for k, bin in self.stream if len(bin) >= self.minAmount}
