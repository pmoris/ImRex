from src.processing.stream import TransformStream


class Labeler(TransformStream):
    def __init__(self, stream, label):
        super().__init__(stream)
        self.label = label

    def transform(self, item, *args, **kwargs):
        return item, self.label


class LabelTrimmer(TransformStream):
    def __init__(self, stream):
        super().__init__(stream)

    def transform(self, item, *args, **kwargs):
        return item[0]
