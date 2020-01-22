import random

from src.bio.peptide_feature import Charge, Hydrophobicity, IsoelectricPoint
from src.bio.util import AMINO_ACIDS


class Peptide(str):
    def __init__(self, seq):
        super().__init__()
        if not (isinstance(seq, Peptide) or isinstance(seq, str)):
            raise RuntimeError("Peptide sequence should be of type str or Peptide")

    def get_charge_vector(self):
        return Charge().calculate(self)

    def get_hydrophobicity_vector(self):
        return Hydrophobicity().calculate(self)

    def get_isoelectric_vector(self):
        return IsoelectricPoint().calculate(self)

    @staticmethod
    def random(length=10):
        return Peptide("".join(random.choice(AMINO_ACIDS) for _ in range(length)))

    def generate_match(self, feature, length=5, random_inverse=True):
        start_index = random.randint(0, len(self) - length)
        match = ""
        for index, aa in enumerate(self):
            if index < start_index:  # start is random
                match += random.choice(AMINO_ACIDS)
            elif index < start_index + length:  # generate random match
                t = feature.generate_match(aa)
                match += t
            else:  # fill with random
                match += random.choice(AMINO_ACIDS)

        if random_inverse and random.randint(1, 2) == 1:
            match = "".join(reversed(match))
        return Peptide(match)
