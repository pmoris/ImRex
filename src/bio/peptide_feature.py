from functools import lru_cache
import random

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData
from Bio.SeqUtils import molecular_weight
import numpy as np
from pyteomics import electrochem

from src.bio.operator import (
    ProductOperator,
    LayeredOperator,
    AbsDifferenceOperator,
    DifferenceOperator,
)

PH = 7


class PeptideFeature(object):
    name = None

    def __repr__(self):
        return self.__class__.__name__ + "()"

    def __init__(self):
        pass

    def _calculate(self, aa):
        """ return the associated value for a given amino acid """
        raise NotImplementedError()

    def generateMatch(self, aa):
        raise RuntimeError("generateMatch not implemented")

    def getBestOperator(self):
        raise NotImplementedError()

    @property
    @lru_cache()
    def values(self):
        return {aa: self._calculate(aa) for aa in AMINO_ACIDS}

    @property
    @lru_cache()
    def max(self):
        return max(self.values.values())

    @property
    @lru_cache()
    def min(self):
        return min(self.values.values())

    def value(self, aa):
        return self.values.get(aa, 0)

    def calculate(self, peptide):
        values = [self.value(aa) for aa in peptide]
        return np.asanyarray(values)

    def matrix(self, pep1, pep2, operator="best"):
        """ By default implements matrix multiplication, but can be overridden """
        if operator == "best":
            operator = self.getBestOperator()
        m = operator.matrix(self.calculate(pep1), self.calculate(pep2))
        return m

    def image_matrix(self, pep1, pep2, operator="best"):
        if operator == "best":
            operator = self.getBestOperator()
        return operator.image_matrix(
            self.calculate(pep1), self.calculate(pep2), self
        ).astype(np.uint8)

    def norm_matrix(self, pep1, pep2, operator="best"):
        if operator == "best":
            operator = self.getBestOperator()
        return operator.norm_matrix(self.calculate(pep1), self.calculate(pep2), self)


@lru_cache()
class Charge(PeptideFeature):
    name = "Charge"

    def __init__(self, ph=PH):
        super().__init__()
        self.ph = ph

    def _calculate(self, aa):
        return electrochem.charge(aa, self.ph)

    def getBestOperator(self):
        return AbsDifferenceOperator()

    # @after(lambda x: print(x, Charge().value(x)))
    def generateMatch(self, amino):
        # print(amino, self.value(amino))

        """ This method uses charge as weight for sampling (neg to pos). """
        # pos, posW = zip(*[(aa, charge) for aa, charge in self.values.items() if charge >= 0])
        # neg, negW = zip(*[(aa, abs(charge)) for aa, charge in self.values.items() if charge < 0])
        # if self.value(aa) >= 0:                  # match pos aa with negative
        #     return random.choices(neg, negW)[0]
        # else:
        #     return random.choices(pos, posW)[0]

        """ This method matches pos to neg, neg to pos and neutral to neutral. """
        CUTOFF = 0.5

        pos = [aa for aa, charge in self.values.items() if charge >= CUTOFF]
        neg = [aa for aa, charge in self.values.items() if charge <= CUTOFF]
        net = [aa for aa, charge in self.values.items() if abs(charge) < CUTOFF]
        if self.value(amino) >= CUTOFF:
            return random.choice(neg)
        elif self.value(amino) <= CUTOFF:
            return random.choice(pos)
        else:
            return random.choice(net)


@lru_cache()
class Hydrophobicity(PeptideFeature):
    """Calculate the hydrophobicity of the amino acids in a protein sequence.

    kd = # Kyte & Doolittle index of hydrophobicity

    Relies on biopython's Bio.SeqUtils.ProtParamData
    Source: https://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParamData-pysrc.html
    """

    name = "Hydrophobicity"

    @property
    @lru_cache()
    def values(self):
        return ProtParamData.kd

    def _calculate(self, aa):
        return ProtParamData.kd[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()

    # @after(lambda x: print(x, Hydrophobicity().value(x), '\n'))
    def generateMatch(self, amino):
        # print(amino, self.value(amino))

        """ This method selects a value close to the value of the aa. """
        acids, weights = zip(
            *[(aa, abs(pow(v - self.value(amino), 2))) for aa, v in self.values.items()]
        )
        secondBest = sorted([w for w in weights if w > 0])[
            0
        ]  # give same amino acid same weight as second best
        invWeights = [1.0 / w if w != 0 else 1 / secondBest for w in weights]
        return random.choices(acids, invWeights)[0]


@lru_cache()
class Polarity(PeptideFeature):
    name = "Polarity"

    def _calculate(self, aa):
        return ProteinAnalysis(aa).isoelectric_point()

    def getBestOperator(self):
        return AbsDifferenceOperator()

    def generateMatch(self, amino):
        # print(amino, self.value(amino))

        """ This method selects a value close to the value of the aa. """
        acids, weights = zip(
            *[(aa, abs(pow(v - self.value(amino), 2))) for aa, v in self.values.items()]
        )
        secondBest = sorted([w for w in weights if w > 0])[
            0
        ]  # give same amino acid same weight as second best
        invWeights = [1.0 / w if w != 0 else 1 / secondBest for w in weights]
        return random.choices(acids, invWeights)[0]


@lru_cache()
class Mass(PeptideFeature):
    name = "Mass"

    def _calculate(self, aa):
        return molecular_weight(
            aa, seq_type="protein", circular=True
        )  # circular to not include water

    def getBestOperator(self):
        return AbsDifferenceOperator()

    def generateMatch(self, amino):
        """ This method selects a value close to the value of the aa. """
        acids, weights = zip(
            *[(aa, abs(pow(v - self.value(amino), 2))) for aa, v in self.values.items()]
        )
        secondBest = sorted([w for w in weights if w > 0])[
            0
        ]  # give same amino acid same weight as second best
        invWeights = [1.0 / w if w != 0 else 1 / secondBest for w in weights]
        return random.choices(acids, invWeights)[0]


@lru_cache()
class Hydrophilicity(PeptideFeature):
    """Calculate the Hydrophilicity of the amino acids in a protein sequence.

    hw = Hopp & Wood Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981). 

    Relies on biopython's Bio.SeqUtils.ProtParamData
    Source: https://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParamData-pysrc.html
    """

    name = "Hydrophilicity"

    @property
    @lru_cache()
    def values(self):
        return ProtParamData.hw

    def _calculate(self, aa):
        return ProtParamData.hw[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class Surface(PeptideFeature):
    """Calculate the surface accessibility of the amino acids in a protein sequence.

    em = Vergoten G & Theophanides T, Biomolecular Structure and Dynamics, pg.138 (1997).

    Relies on biopython's Bio.SeqUtils.ProtParamData
    Source: https://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParamData-pysrc.html
    """

    name = "SurfaceAccessibility"

    @property
    @lru_cache()
    def values(self):
        return ProtParamData.em

    def _calculate(self, aa):
        return ProtParamData.em[aa]


@lru_cache()
class Flexibility(PeptideFeature):
    """Calculate the flexibility of the amino acids in a protein sequence.

    Flex = Normalized flexibility parameters (B-values), average Vihinen M., Torkkila E., Riikonen P. Proteins. 19(2):141-9(1994).

    Relies on biopython's Bio.SeqUtils.ProtParamData
    Source: https://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParamData-pysrc.html
    """

    name = "Flexibility"

    @property
    @lru_cache()
    def values(self):
        return ProtParamData.Flex

    def _calculate(self, aa):
        return ProtParamData.Flex[aa]


@lru_cache()
class Transfer(PeptideFeature):
    """Calculate the surface transfer energy of the amino acids in a protein sequence.

    ja = 2 Janin Interior to surface transfer energy scale

    Relies on biopython's Bio.SeqUtils.ProtParamData
    Source: https://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParamData-pysrc.html
    """

    name = "InteriorToSurfaceTransferEnergy"

    @property
    @lru_cache()
    def values(self):
        return ProtParamData.ja

    def _calculate(self, aa):
        return ProtParamData.ja[aa]


@lru_cache()
class Prime(PeptideFeature):
    name = "Prime"

    @property
    @lru_cache()
    def values(self):
        return {aa: float(p) for aa, p in zip(AMINO_ACIDS, gen_primes())}

    def getBestOperator(self):
        return ProductOperator()


# Source: https://github.com/bittremieux/TCR-Classifier/blob/master/tcr_classifier_v2.ipynb
basicity = {
    "A": 206.4,
    "B": 210.7,
    "C": 206.2,
    "D": 208.6,
    "E": 215.6,
    "F": 212.1,
    "G": 202.7,
    "H": 223.7,
    "I": 210.8,
    "K": 221.8,
    "L": 209.6,
    "M": 213.3,
    "N": 212.8,
    "P": 214.4,
    "Q": 214.2,
    "R": 237.0,
    "S": 207.6,
    "T": 211.7,
    "V": 208.7,
    "W": 216.1,
    "X": 210.2,
    "Y": 213.1,
    "Z": 214.9,
}

hydrophobicity = {
    "A": 0.16,
    "B": -3.14,
    "C": 2.50,
    "D": -2.49,
    "E": -1.50,
    "F": 5.00,
    "G": -3.31,
    "H": -4.63,
    "I": 4.41,
    "K": -5.00,
    "L": 4.76,
    "M": 3.23,
    "N": -3.79,
    "P": -4.92,
    "Q": -2.76,
    "R": -2.77,
    "S": -2.85,
    "T": -1.08,
    "V": 3.02,
    "W": 4.88,
    "X": 4.59,
    "Y": 2.00,
    "Z": -2.13,
}

helicity = {
    "A": 1.24,
    "B": 0.92,
    "C": 0.79,
    "D": 0.89,
    "E": 0.85,
    "F": 1.26,
    "G": 1.15,
    "H": 0.97,
    "I": 1.29,
    "K": 0.88,
    "L": 1.28,
    "M": 1.22,
    "N": 0.94,
    "P": 0.57,
    "Q": 0.96,
    "R": 0.95,
    "S": 1.00,
    "T": 1.09,
    "V": 1.27,
    "W": 1.07,
    "X": 1.29,
    "Y": 1.11,
    "Z": 0.91,
}

mutation_stability = {
    "A": 13.0,
    "C": 52.0,
    "D": 11.0,
    "E": 12.0,
    "F": 32.0,
    "G": 27.0,
    "H": 15.0,
    "I": 10.0,
    "K": 24.0,
    "L": 34.0,
    "M": 6.0,
    "N": 6.0,
    "P": 20.0,
    "Q": 10.0,
    "R": 17.0,
    "S": 10.0,
    "T": 11.0,
    "V": 17.0,
    "W": 55.0,
    "Y": 31.0,
}


@lru_cache()
class TCRexBasicity(PeptideFeature):
    name = "TCRexBasicity"

    @property
    @lru_cache()
    def values(self):
        return basicity

    def _calculate(self, aa):
        return basicity[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class TCRexHydrophobicity(PeptideFeature):
    name = "TCRexHydrophobicity"

    @property
    @lru_cache()
    def values(self):
        return hydrophobicity

    def _calculate(self, aa):
        return hydrophobicity[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class TCRexHelicity(PeptideFeature):
    name = "TCRexHelicity"

    @property
    @lru_cache()
    def values(self):
        return helicity

    def _calculate(self, aa):
        return helicity[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class TCRexMutationStability(PeptideFeature):
    name = "TCRexMutationStability"

    @property
    @lru_cache()
    def values(self):
        return mutation_stability

    def _calculate(self, aa):
        return mutation_stability[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


# Sieve of Eratosthenes
# Code by David Eppstein, UC Irvine, 28 Feb 2002
# http://code.activestate.com/recipes/117119/
def gen_primes():
    """ Generate an infinite sequence of prime numbers. """
    # Maps composites to primes witnessing their compositeness.
    # This is memory efficient, as the sieve is not "run forward"
    # indefinitely, but only as long as required by the current
    # number being tested.
    #
    D = {}

    # The running integer that's checked for primeness
    q = 2

    while True:
        if q not in D:
            # q is a new prime.
            # Yield it and mark its first multiple that isn't
            # already marked in previous iterations
            #
            yield q
            D[q * q] = [q]
        else:
            # q is composite. D[q] is the list of primes that
            # divide it. Since we've reached q, we no longer
            # need it in the map, but we'll mark the next
            # multiples of its witnesses to prepare for larger
            # numbers
            #
            for p in D[q]:
                D.setdefault(p + q, []).append(p)
            del D[q]

        q += 1


featuresMap = {
    "charge": Charge(),
    "hydrophob": Hydrophobicity(),
    "hydrophil": Hydrophilicity(),
    "polarity": Polarity(),
    "mass": Mass(),
    "prime": Prime(),
    "basicity": TCRexBasicity(),
    "helicity": TCRexHelicity(),
    "hydrophob2": TCRexHydrophobicity(),
    "mutationstab": TCRexMutationStability(),
}

operatorsMap = {
    "prod": ProductOperator(),
    "absdiff": AbsDifferenceOperator(),
    "diff": DifferenceOperator(),
    "layer": LayeredOperator(),
    "best": "best",
}


def parseFeatures(string):
    names = [name.strip() for name in string.split(",")]
    try:
        return [featuresMap[name] for name in names]
    except ValueError as e:
        print("Unkown feature encountered")
        raise e


def parseOperator(string):
    return operatorsMap[string]
