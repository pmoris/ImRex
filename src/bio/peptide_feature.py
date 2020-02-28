from functools import lru_cache
import random

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils import ProtParamData
import numpy as np
from pyteomics import electrochem

from src.bio.operator import (
    ProductOperator,
    LayeredOperator,
    AbsDifferenceOperator,
    DifferenceOperator,
)
from src.definitions.amino_acid_properties import (
    AMINO_ACIDS,
    ATCHLEY_FACTOR_1,
    ATCHLEY_FACTOR_2,
    ATCHLEY_FACTOR_3,
    ATCHLEY_FACTOR_4,
    ATCHLEY_FACTOR_5,
    TCREX_BASICITY,
    TCREX_HYDROPHOBICITY,
    TCREX_HELICITY,
    TCREX_MUTATION_STABILITY,
    PH,
)


class PeptideFeature(object):
    name = None

    def __repr__(self):
        return self.__class__.__name__ + "()"

    def __init__(self):
        pass

    def _calculate(self, aa):
        """ Return the associated value of the property for a given amino acid.

        Implemented by the feature subclasses.

        Parameters
        ----------
        aa : string
            amino acid

        Returns
        -------
        float or int

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError()

    def generateMatch(self, aa):
        raise RuntimeError("generateMatch not implemented")

    def getBestOperator(self):
        raise NotImplementedError()

    def get_max_feature_value(self):
        raise NotImplementedError()

    def get_min_feature_value(self):
        raise NotImplementedError()

    @property
    @lru_cache()
    def values(self):
        """Return a dictionary mapping all amino acids to the property.

        Most feature subclasses override this method and thus skip the _calculate() call.

        Returns
        -------
        dict
        """
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
        """Return the associated value of the property for all amino acids in a given peptide sequence,
        and 0 when the amino acid is not found or invalid.

        Is never overriden by a feature subclass.

        Parameters
        ----------
        peptide : string
            A peptide sequence for which the property should be calculated per amino acid.

        Returns
        -------
        array
            An array of property values (float/int).
        """
        # values = [self.value(aa) for aa in peptide]
        values = [self.values.get(aa, 0) for aa in peptide]
        # values = [self._calculate(aa) for aa in peptide]
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
    # name = "Polarity"
    name = "Isoelectric point"

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


@lru_cache()
class TCRexBasicity(PeptideFeature):
    name = "TCRexBasicity"

    @property
    @lru_cache()
    def values(self):
        return TCREX_BASICITY

    def _calculate(self, aa):
        return TCREX_BASICITY[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class TCRexHydrophobicity(PeptideFeature):
    name = "TCRexHydrophobicity"

    @property
    @lru_cache()
    def values(self):
        return TCREX_HYDROPHOBICITY

    def _calculate(self, aa):
        return TCREX_HYDROPHOBICITY[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class TCRexHelicity(PeptideFeature):
    name = "TCRexHelicity"

    @property
    @lru_cache()
    def values(self):
        return TCREX_HELICITY

    def _calculate(self, aa):
        return TCREX_HELICITY[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class TCRexMutationStability(PeptideFeature):
    name = "TCRexMutationStability"

    @property
    @lru_cache()
    def values(self):
        return TCREX_MUTATION_STABILITY

    def _calculate(self, aa):
        return TCREX_MUTATION_STABILITY[aa]

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


@lru_cache()
class AtchleyFactor1(PeptideFeature):
    name = "Atchley_factor_1"

    @property
    @lru_cache()
    def values(self):
        return ATCHLEY_FACTOR_1

    def _calculate(self, aa):
        return ATCHLEY_FACTOR_1[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class AtchleyFactor2(PeptideFeature):
    name = "Atchley_factor_2"

    @property
    @lru_cache()
    def values(self):
        return ATCHLEY_FACTOR_2

    def _calculate(self, aa):
        return ATCHLEY_FACTOR_2[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class AtchleyFactor3(PeptideFeature):
    name = "Atchley_factor_3"

    @property
    @lru_cache()
    def values(self):
        return ATCHLEY_FACTOR_3

    def _calculate(self, aa):
        return ATCHLEY_FACTOR_3[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class AtchleyFactor4(PeptideFeature):
    name = "Atchley_factor_4"

    @property
    @lru_cache()
    def values(self):
        return ATCHLEY_FACTOR_4

    def _calculate(self, aa):
        return ATCHLEY_FACTOR_4[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


@lru_cache()
class AtchleyFactor5(PeptideFeature):
    name = "Atchley_factor_5"

    @property
    @lru_cache()
    def values(self):
        return ATCHLEY_FACTOR_5

    def _calculate(self, aa):
        return ATCHLEY_FACTOR_5[aa]

    def getBestOperator(self):
        return AbsDifferenceOperator()


featuresMap = {
    "charge": Charge(),
    "hydrophob": Hydrophobicity(),
    "hydrophil": Hydrophilicity(),
    "polarity": Polarity(),
    "mass": Mass(),
    "surface": Surface(),
    "flexibility": Flexibility(),
    "transfer": Transfer(),
    "prime": Prime(),
    "basicity": TCRexBasicity(),
    "helicity": TCRexHelicity(),
    "hydrophob2": TCRexHydrophobicity(),
    "mutationstab": TCRexMutationStability(),
    "atchley1": AtchleyFactor1(),
    "atchley2": AtchleyFactor2(),
    "atchley3": AtchleyFactor3(),
    "atchley4": AtchleyFactor4(),
    "atchley5": AtchleyFactor5(),
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
