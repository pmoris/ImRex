import numpy as np
from PIL import Image


def image_from_matrix(matrix, mode="L", index=None):
    """ Generate image from matrix with values between 0 and 255. """
    if mode == "L":
        return Image.fromarray(matrix.astype(np.uint8), "L")
    elif mode in ["RGB", "CMYK"]:
        amt_channels = len(mode)
        img = np.zeros_like(matrix, dtype=np.uint8)
        img = np.expand_dims(img, axis=-1)
        img = np.repeat(img, amt_channels, axis=2)
        img[..., index] = matrix
        return Image.fromarray(img, mode)
    else:
        raise RuntimeError("Unkown image mode")


def image_from_matrices(*matrices, mode="RGB"):
    """ Generate rgb image from matrices with values between 0 and 255. """
    return image_from_tensor(np.dstack(matrices).astype(np.uint8), mode=mode)


def image_from_tensor(tensor, mode="RGB"):
    """ Generate rgb image from tensor of 3 matrices with values between 0 and 255. """
    img = Image.fromarray(tensor.astype(np.uint8), mode=mode)
    return img
