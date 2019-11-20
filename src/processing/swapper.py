import random

from src.processing.stream import TransformStream


class Swapper(TransformStream):
    """ Swap two inputs. """

    def __init__(self, stream):
        super().__init__(stream)

    def transform(self, item, *args, **kwargs):
        x, y = item
        if random.randint(0, 1) == 0:
            return x, y
        else:
            return reversed(x), y
