from src.processing.stream import TransformStream


class ImageGenerator(TransformStream):
    def __init__(self, stream, feature_builder, has_label=True):
        super().__init__(stream)
        self.feature_builder = feature_builder
        self.has_label = has_label

    def transform(self, item, *args, **kwargs):
        if self.has_label:
            x, y = item
            return self.feature_builder.generate_feature(x), y
        else:
            return self.feature_builder.generate_feature(item)
