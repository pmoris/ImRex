from src.processing.data_stream import DataStream
from src.processing.stream import Stream


class Filter(Stream):
    def __init__(self, stream: DataStream, has_label=True):
        """Parent class for filtering DataStreams based on a certain property.

        Note that when this Filter object is iterated over, the original DataStream
        will be consumed. This will cause get() or next() calls to
        raise a StopIteration error, but will not prevent loops or lists
        from being made.

        Note that the Filter object itself will consequently be completely
        consumed as well, and even loops or list calls will fail. Adding
        an __iter__ call to self.stream would prevent this.

        Parameters
        ----------
        stream : DataStream
            The DataStream to be filtered.
        has_label : bool, optional
            Whether each element of the DataStream contains a class label or not, by default True.
            E.g. (('CASSFQGSNEKLFF', 'RAKFKQLL'), 1) versus 'CASSFQGSNEKLFF'
        """
        super().__init__(stream)
        self.stream = stream
        self.has_label = has_label

    def __len__(self):
        return None  # size not needed

    # def __iter__(self):
    #     return iter(self.stream)

    def is_valid(self, item):
        raise NotImplementedError

    def get(self, *args, **kwargs):
        while True:
            # get next element from stream object
            element = self.stream.get(*args, **kwargs)
            # if the stream object has a label, retrieve the sequence pair
            if self.has_label:
                item = element[0]
            # else, just use the element as a whole
            else:
                item = element
            # return the element if it satisfies the criteria, e.g. length
            if self.is_valid(item):
                return element


class SizeFilter(Filter):
    def __init__(self, stream: DataStream, *ranges, has_label=True):
        """ Filter DataStreams based on the length of its elements.

        When applied to cdr3 reference sequence DataStreams, has_label should be set to False.

        Parameters
        ----------
        ranges : Tuple
            The range restrictions for the elements in the DataStream.
            Format: ((10, 20), (8, 13))
            Where the first tuple applies to the first (cdr3) sequence,
            and the second to the second (epitope) sequence.
        """
        super().__init__(stream, has_label=has_label)
        self.ranges = ranges

    def is_valid(self, item):
        # if the stream contains a single element per observation e.g. a cdr3 reference sequence
        if len(self.ranges) == 1 and not isinstance(item, list):
            # convert it into the first element of a tuple
            # this will ensure that the first element is checked against
            # the first length range pair, and the second length range pair is ignored
            item = (item,)

        # for each sequence element of the tuple (or sequence and None in case of cd3r refs)
        # check if it falls within the length ranges
        # ranges contains the cdr3 ranges first and the epitope ranges second
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
