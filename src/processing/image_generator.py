from src.processing.stream import TransformStream


class ImageGenerator(TransformStream):
    def __init__(self, stream, featureBuilder):
        super().__init__(stream)
        self.featureBuilder = featureBuilder

    def transform(self, item, *args, **kwargs):
        x, y = item
        return self.featureBuilder.generateFeature(x), y
