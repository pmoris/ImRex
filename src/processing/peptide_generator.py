import random

from .stream import Stream
from bio.peptide import Peptide


class PeptideGenerator(Stream):
    def __init__(self, minSize, maxSize):
        super().__init__()
        self.minSize = minSize
        self.maxSize = maxSize

    def __len__(self):
        return None

    def get(self, *args, **kwargs):
        length = random.randint(self.minSize, self.maxSize)
        return Peptide.random(length)
