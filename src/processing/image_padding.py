import numpy as np

from src.processing.stream import TransformStream


class ImagePadding(TransformStream):
    def __init__(self, stream, width, height, pad_value=0, has_label=True):
        super().__init__(stream)
        self.width = width
        self.height = height
        self.pad_value = pad_value
        self.has_label = has_label

    def transform(self, item, *args, **kwargs):
        if self.has_label:
            image, label = item
            return self._padding(image), label
        else:
            return self._padding(item)

    def _padding(self, image):
        hor_padding = self.width - image.shape[0]
        ver_padding = self.height - image.shape[1]

        hor_padding_before = int(hor_padding // 2)
        hor_padding_after = hor_padding - hor_padding_before

        ver_padding_before = int(ver_padding // 2)
        ver_padding_after = ver_padding - ver_padding_before

        padding = (
            (hor_padding_before, hor_padding_after),
            (ver_padding_before, ver_padding_after),
            (0, 0),
        )

        padded = np.pad(image, padding, mode="constant", constant_values=self.pad_value)

        return padded
