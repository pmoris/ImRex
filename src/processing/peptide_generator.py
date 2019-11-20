import random

from src.processing.stream import Stream
from src.bio.peptide import Peptide


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
