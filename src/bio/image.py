import numpy as np
from PIL import Image


def imageFromMatrix(matrix, mode="L", index=None):
    """ Generate image from matrix with values between 0 and 255 """
    if mode == "L":
        return Image.fromarray(matrix.astype(np.uint8), "L")
    elif mode in ["RGB", "CMYK"]:
        amtChannels = len(mode)
        img = np.zeros_like(matrix, dtype=np.uint8)
        img = np.expand_dims(img, axis=-1)
        img = np.repeat(img, amtChannels, axis=2)
        img[..., index] = matrix
        return Image.fromarray(img, mode)
    else:
        raise RuntimeError("Unkown image mode")


def imageFromMatrices(*matrices, mode="RGB"):
    """ Generate rgb image from matrices with values between 0 and 255 """
    return imageFromTensor(np.dstack(matrices).astype(np.uint8), mode=mode)


def imageFromTensor(tensor, mode="RGB"):
    """ Generate rgb image from tensor of 3 matrices with values between 0 and 255 """
    img = Image.fromarray(tensor.astype(np.uint8), mode=mode)
    return img
