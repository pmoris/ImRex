import numpy as np

from src.processing.stream import TransformStream


class ImagePadding(TransformStream):
    def __init__(self, stream, width, height, padValue=0):
        super().__init__(stream)
        self.width = width
        self.height = height
        self.padValue = padValue

    def transform(self, item, *args, **kwargs):
        image, label = item
        horPadding = self.width - image.shape[0]
        verPadding = self.height - image.shape[1]

        horPadding1 = int(horPadding // 2)
        horPadding2 = horPadding - horPadding1

        verPadding1 = int(verPadding // 2)
        verPadding2 = verPadding - verPadding1

        padding = ((horPadding1, horPadding2), (verPadding1, verPadding2), (0, 0))
        # print(image)
        # print(image.shape)
        # print(padding)
        padded = np.pad(image, padding, mode="constant", constant_values=self.padValue)
        # print("Padded an image")
        return padded, label
