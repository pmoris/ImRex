from typing import Optional

from src.processing.stream import TransformStream


class Reverser(TransformStream):
    """ Reverse input sequences. """

    def __init__(self, stream, target: Optional[str] = None):
        """ Reverses the sequences in the DataStream.

        Parameters
        ----------
        stream : DataStream
            The DataStream whose sequences should be reversed.
        target : str
            Indicates which sequence should be reversed. Either "cdr3" or "epitope".
            If omitted reversed both sequences.
        """
        super().__init__(stream)
        self.target = target

    def transform(self, item, *args, **kwargs):
        sequences, label = item
        cdr3, epitope = sequences

        if self.target == "cdr3":
            cdr3 = cdr3[::-1]
        elif self.target == "epitope":
            epitope = epitope[::-1]
        else:
            cdr3 = cdr3[::-1]
            epitope = epitope[::-1]

        return ((cdr3, epitope), label)
