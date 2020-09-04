from src.processing.stream import TransformStream


class Decoyer(TransformStream):
    """ Transform second (=epitope) sequence into decoy. """

    def __init__(self, stream, decoy_dict: dict):
        super().__init__(stream)
        self.decoy_dict = decoy_dict

    def transform(self, item, *args, **kwargs):
        sequence_tuple, label = item
        cdr3, epitope = sequence_tuple
        return [(cdr3, self.decoy_dict[epitope]), label]
