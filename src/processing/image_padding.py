import numpy as np

from src.processing.stream import TransformStream


class ImagePadding(TransformStream):
    def __init__(self, stream, width, height, pad_value=0):
        super().__init__(stream)
        self.width = width
        self.height = height
        self.pad_value = pad_value

    def transform(self, item, *args, **kwargs):
        image, label = item
        hor_padding = self.width - image.shape[0]
        ver_padding = self.height - image.shape[1]

        hor_padding1 = int(hor_padding // 2)
        hor_padding2 = hor_padding - hor_padding1

        ver_padding1 = int(ver_padding // 2)
        ver_padding2 = ver_padding - ver_padding1

        padding = ((hor_padding1, hor_padding2), (ver_padding1, ver_padding2), (0, 0))
        # print(image)
        # print(image.shape)
        # print(padding)
        padded = np.pad(image, padding, mode="constant", constant_values=self.pad_value)
        # print("Padded an image")
        return padded, label
