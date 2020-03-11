import random

from src.bio.peptide import Peptide
from src.processing.stream import Stream


class PeptideGenerator(Stream):
    def __init__(self, min_size, max_size):
        super().__init__()
        self.min_size = min_size
        self.max_size = max_size

    def __len__(self):
        return None

    def get(self, *args, **kwargs):
        length = random.randint(self.min_size, self.max_size)
        return Peptide.random(length)
