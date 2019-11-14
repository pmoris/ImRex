import random

from bio.peptide_feature import Charge, Hydrophobicity, Polarity
from bio.util import AMINO_ACIDS


class Peptide(str):
    def __init__(self, seq):
        super().__init__()
        if not (isinstance(seq, Peptide) or isinstance(seq, str)):
            raise RuntimeError("Peptide sequence should be of type str or Peptide")

    def getChargeVector(self):
        return Charge().calculate(self)

    def getHyrdophobicityVector(self):
        return Hydrophobicity().calculate(self)

    def getPolarityVector(self):
        return Polarity().calculate(self)

    @staticmethod
    def random(length=10):
        return Peptide(''.join(random.choice(AMINO_ACIDS) for _ in range(length)))

    def generateMatch(self, feature, length=5, random_inverse=True):
        startIndex = random.randint(0, len(self) - length)
        match = ""
        for index, aa in enumerate(self):
            if index < startIndex:                      # start is random
                match += random.choice(AMINO_ACIDS)
            elif index < startIndex + length:           # generate random match
                t = feature.generateMatch(aa)
                match += t
            else:                                       # fill with random
                match += random.choice(AMINO_ACIDS)

        if random_inverse and random.randint(1, 2) == 1:
            match = ''.join(reversed(match))
        return Peptide(match)

