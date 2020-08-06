from collections import defaultdict
from functools import lru_cache

from src.processing.stream import GroupedStream


class Grouper(GroupedStream):
    def __init__(self, stream):
        super().__init__(stream)
        self.stream = stream

    def __len__(self):
        return len(self.stream)

    def get_key(self, item):
        raise NotImplementedError

    @lru_cache()
    def get_groups(self):
        bins = defaultdict(list)
        for item in self.stream:
            key = self.get_key(item)
            bins[key].append(item)
        return bins


class ShapeGrouper(Grouper):
    def __init__(self, stream, contains_label=True):
        super().__init__(stream)
        self.contains_label = contains_label

    def get_key(self, item):
        if self.contains_label:
            item = item[0]
        return tuple(len(i) for i in item)


class SizeGrouper(Grouper):
    def __init__(self, stream, contains_label=True):
        super().__init__(stream)
        self.contains_label = contains_label

    def get_key(self, item):
        if self.contains_label:
            item = item[0]
        return len(item)


class GroupedAmountFilter(GroupedStream):
    def __init__(self, stream, min_amount):
        super().__init__(stream)
        self.stream = stream
        self.min_amount = min_amount

    def __len__(self):
        return sum([len(binned) for binned in self.get_groups().values()])

    @lru_cache()
    def get_groups(self):
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

        return {k: bin for k, bin in self.stream if len(bin) >= self.min_amount}
