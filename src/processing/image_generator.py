from src.processing.stream import TransformStream


class ImageGenerator(TransformStream):
    def __init__(self, stream, feature_builder):
        super().__init__(stream)
        self.feature_builder = feature_builder

    def transform(self, item, *args, **kwargs):
        x, y = item
        return self.feature_builder.generate_feature(x), y
