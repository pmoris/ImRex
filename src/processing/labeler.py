from src.processing.data_stream import DataStream
from src.processing.stream import TransformStream


class Labeler(TransformStream):
    def __init__(self, stream, label):
        super().__init__(stream)
        self.label = label

    def transform(self, item, *args, **kwargs):
        return item, self.label


class LabelTrimmer(TransformStream):
    """Stream object that removes the class label from the input.

    More specifically, it transforms the provided Stream by keeping
    only the first element of each observation.

    """

    def __init__(self, stream: DataStream):
        super().__init__(stream)

    def transform(self, item, *args, **kwargs):
        """Trims the class label from each observation.

        Parameters
        ----------
        item : Tuple
            (('CAVLSLSGSARQLTF', 'NLVPMVATV'), 1)

        Returns
        -------
        Tuple
            ('CAVLSLSGSARQLTF', 'NLVPMVATV')
        """
        return item[0]
