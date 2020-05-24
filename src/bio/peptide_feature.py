from functools import lru_cache
import random

from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils import ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
from pyteomics import electrochem

from src.bio.operator import (
    AbsDifferenceOperator,
    DifferenceOperator,
    LayeredOperator,
    ProductOperator,
)
from src.definitions.amino_acid_properties import (
    AMINO_ACIDS,
    ATCHLEY_FACTOR_1,
    ATCHLEY_FACTOR_2,
    ATCHLEY_FACTOR_3,
    ATCHLEY_FACTOR_4,
    ATCHLEY_FACTOR_5,
    PH,
    TCREX_BASICITY,
    TCREX_HELICITY,
    TCREX_HYDROPHOBICITY,
    TCREX_MUTATION_STABILITY,
    KIDERA_FACTOR_1,
    KIDERA_FACTOR_2,
    KIDERA_FACTOR_3,
    KIDERA_FACTOR_4,
    KIDERA_FACTOR_5,
    KIDERA_FACTOR_6,
    KIDERA_FACTOR_7,
    KIDERA_FACTOR_8,
    KIDERA_FACTOR_9,
    KIDERA_FACTOR_10,
)


class PeptideFeature(object):
    name = None

    def __repr__(self):
        return self.__class__.__name__ + "()"

    def __init__(self):
        pass

    def _calculate(self, aa: str) -> float:
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

    def generate_match(self, aa):
        raise RuntimeError("generateMatch not implemented")

    def get_best_operator(self):
        return AbsDifferenceOperator()

    def get_max_feature_value(self):
        raise NotImplementedError()

    def get_min_feature_value(self):
        raise NotImplementedError()

    @property
    @lru_cache()
    def values(self) -> dict:
        """Return a dictionary mapping all amino acids to the property.

        Most feature subclasses override this method and thus skip the _calculate() call.

        Returns
        -------
        dict
        """
        return {aa: self._calculate(aa) for aa in AMINO_ACIDS}

    @property  # noqa: A003
    @lru_cache()
    def max(self):
        return max(self.values.values())

    @property  # noqa: A003
    @lru_cache()
    def min(self):
        return min(self.values.values())

    def calculate(self, peptide: str) -> np.ndarray:
        """Return the associated value of the property for all amino acids in a given peptide sequence, and 0 when the amino acid is not found or invalid.

        Is never overriden by a feature subclass.

        Parameters
        ----------
        peptide : string
            A peptide sequence for which the property should be calculated per amino acid.

        Returns
        -------
        array
            An array of property values (type = np.float64).
        """
        values = [self.values.get(aa, 0) for aa in peptide]
        # or equivalently: values = [self._calculate(aa) for aa in peptide]
        return np.asanyarray(values)

    def matrix(self, pep1: str, pep2: str, operator="best") -> np.ndarray:
        """Compute the pairwise amino acid matrix of two amino acid sequences using the given operator, for this amino acid property/feature.

        Parameters
        ----------
        pep1 : str
            First amino acid sequence
        pep2 : str
            Second amino acid sequence
        operator : str, optional
            The operator to use for combining the amino acid properties, by default "best".

        Returns
        -------
        np.ndarray
            The combined pairwise amino acid properties in matrix format.
        """
        if operator == "best":
            operator = self.get_best_operator()
        m = operator.matrix(self.calculate(pep1), self.calculate(pep2))
        return m

    def image_matrix(self, pep1: str, pep2: str, operator="best") -> np.ndarray:
        """Compute the scaled pairwise amino acid matrix of two amino acid sequences using the given operator, for this amino acid property/feature.

        Elements are scaled between 0 and 1, where the minimum and maximum value
        are defined by the smallest and largest value that the pairwise combination
        of the amino acid property can take among all 20 amino acids
        (i.e. not the minimum and maximum in the matrix under consideration).

        Parameters
        ----------
        pep1 : str
            First amino acid sequence
        pep2 : str
            Second amino acid sequence
        operator : str, optional
            The operator to use for combining the amino acid properties, by default "best".

        Returns
        -------
        np.ndarray
            The combined scaled (0-255) pairwise amino acid properties in matrix format.
        """
        if operator == "best":
            operator = self.get_best_operator()
        return operator.image_matrix(
            self.calculate(pep1), self.calculate(pep2), self
        ).astype(np.uint8)

    def norm_matrix(self, pep1: str, pep2: str, operator="best") -> np.ndarray:
        """Compute the normalized pairwise amino acid matrix of two amino acid sequences using the given operator, for this amino acid property/feature.

        Elements are normalized between 0 and 1, where the minimum and maximum value
        are defined by the smallest and largest value that the pairwise combination
        of the amino acid property can take among all 20 amino acids
        (i.e. not the minimum and maximum in the matrix under consideration).

        Parameters
        ----------
        pep1 : str
            First amino acid sequence
        pep2 : str
            Second amino acid sequence
        operator : str, optional
            The operator to use for combining the amino acid properties, by default "best".

        Returns
        -------
        np.ndarray
            The combined normalized (0-1) pairwise amino acid properties in matrix format.
        """
        if operator == "best":
            operator = self.get_best_operator()
        return operator.norm_matrix(self.calculate(pep1), self.calculate(pep2), self)


@lru_cache()
class Charge(PeptideFeature):
    name = "Charge"

    def __init__(self, ph=PH):
        super().__init__()
        self.ph = ph

    def _calculate(self, aa: str) -> float:
        return electrochem.charge(aa, self.ph)

    def generate_match(self, amino):
        """ Match pos to neg, neg to pos and neutral to neutral. """
        CUTOFF = 0.5

        pos = [aa for aa, charge in self.values.items() if charge >= CUTOFF]
        neg = [aa for aa, charge in self.values.items() if charge <= CUTOFF]
        net = [aa for aa, charge in self.values.items() if abs(charge) < CUTOFF]
        if self._calculate(amino) >= CUTOFF:
            return random.choice(neg)
        elif self._calculate(amino) <= CUTOFF:
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
    def values(self) -> dict:
        return ProtParamData.kd

    def _calculate(self, aa: str) -> float:
        return ProtParamData.kd[aa]

    def generate_match(self, amino):
        """ Select a value close to the value of the aa. """
        acids, weights = zip(
            *[
                (aa, abs(pow(v - self._calculate(amino), 2)))
                for aa, v in self.values.items()
            ]
        )
        second_best = sorted([w for w in weights if w > 0])[
            0
        ]  # give same amino acid same weight as second best
        inv_weights = [1.0 / w if w != 0 else 1 / second_best for w in weights]
        return random.choices(acids, inv_weights)[0]


@lru_cache()
class IsoelectricPoint(PeptideFeature):
    name = "Isoelectric point"

    def _calculate(self, aa: str) -> float:
        return ProteinAnalysis(aa).isoelectric_point()

    def generate_match(self, amino):
        """ Select a value close to the value of the aa. """
        acids, weights = zip(
            *[
                (aa, abs(pow(v - self._calculate(amino), 2)))
                for aa, v in self.values.items()
            ]
        )
        second_best = sorted([w for w in weights if w > 0])[
            0
        ]  # give same amino acid same weight as second best
        inv_weights = [1.0 / w if w != 0 else 1 / second_best for w in weights]
        return random.choices(acids, inv_weights)[0]


@lru_cache()
class Mass(PeptideFeature):
    name = "Mass"

    def _calculate(self, aa: str) -> float:
        return molecular_weight(
            aa, seq_type="protein", circular=True
        )  # circular to not include water

    def generate_match(self, amino):
        """ Select a value close to the value of the aa. """
        acids, weights = zip(
            *[
                (aa, abs(pow(v - self._calculate(amino), 2)))
                for aa, v in self.values.items()
            ]
        )
        second_best = sorted([w for w in weights if w > 0])[
            0
        ]  # give same amino acid same weight as second best
        inv_weights = [1.0 / w if w != 0 else 1 / second_best for w in weights]
        return random.choices(acids, inv_weights)[0]


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
    def values(self) -> dict:
        return ProtParamData.hw

    def _calculate(self, aa: str) -> float:
        return ProtParamData.hw[aa]


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
    def values(self) -> dict:
        return ProtParamData.em

    def _calculate(self, aa: str) -> float:
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
    def values(self) -> dict:
        return ProtParamData.Flex

    def _calculate(self, aa: str) -> float:
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
    def values(self) -> dict:
        return ProtParamData.ja

    def _calculate(self, aa: str) -> float:
        return ProtParamData.ja[aa]


@lru_cache()
class Prime(PeptideFeature):
    name = "Prime"

    @property
    @lru_cache()
    def values(self) -> dict:
        return {aa: float(p) for aa, p in zip(AMINO_ACIDS, gen_primes())}

    def get_best_operator(self):
        return ProductOperator()


@lru_cache()
class TCRexBasicity(PeptideFeature):
    name = "TCRexBasicity"

    @property
    @lru_cache()
    def values(self) -> dict:
        return TCREX_BASICITY

    def _calculate(self, aa: str) -> float:
        return TCREX_BASICITY[aa]


@lru_cache()
class TCRexHydrophobicity(PeptideFeature):
    name = "TCRexHydrophobicity"

    @property
    @lru_cache()
    def values(self) -> dict:
        return TCREX_HYDROPHOBICITY

    def _calculate(self, aa: str) -> float:
        return TCREX_HYDROPHOBICITY[aa]


@lru_cache()
class TCRexHelicity(PeptideFeature):
    name = "TCRexHelicity"

    @property
    @lru_cache()
    def values(self) -> dict:
        return TCREX_HELICITY

    def _calculate(self, aa: str) -> float:
        return TCREX_HELICITY[aa]


@lru_cache()
class TCRexMutationStability(PeptideFeature):
    name = "TCRexMutationStability"

    @property
    @lru_cache()
    def values(self) -> dict:
        return TCREX_MUTATION_STABILITY

    def _calculate(self, aa: str) -> float:
        return TCREX_MUTATION_STABILITY[aa]


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
    def values(self) -> dict:
        return ATCHLEY_FACTOR_1

    def _calculate(self, aa: str) -> float:
        return ATCHLEY_FACTOR_1[aa]


@lru_cache()
class AtchleyFactor2(PeptideFeature):
    name = "Atchley_factor_2"

    @property
    @lru_cache()
    def values(self) -> dict:
        return ATCHLEY_FACTOR_2

    def _calculate(self, aa: str) -> float:
        return ATCHLEY_FACTOR_2[aa]


@lru_cache()
class AtchleyFactor3(PeptideFeature):
    name = "Atchley_factor_3"

    @property
    @lru_cache()
    def values(self) -> dict:
        return ATCHLEY_FACTOR_3

    def _calculate(self, aa: str) -> float:
        return ATCHLEY_FACTOR_3[aa]


@lru_cache()
class AtchleyFactor4(PeptideFeature):
    name = "Atchley_factor_4"

    @property
    @lru_cache()
    def values(self) -> dict:
        return ATCHLEY_FACTOR_4

    def _calculate(self, aa: str) -> float:
        return ATCHLEY_FACTOR_4[aa]


@lru_cache()
class AtchleyFactor5(PeptideFeature):
    name = "Atchley_factor_5"

    @property
    @lru_cache()
    def values(self) -> dict:
        return ATCHLEY_FACTOR_5

    def _calculate(self, aa: str) -> float:
        return ATCHLEY_FACTOR_5[aa]


@lru_cache()
class KideraFactor1(PeptideFeature):
    name = "Kidera_factor_1"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_1

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_1[aa]


@lru_cache()
class KideraFactor2(PeptideFeature):
    name = "Kidera_factor_2"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_2

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_2[aa]


@lru_cache()
class KideraFactor3(PeptideFeature):
    name = "Kidera_factor_3"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_3

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_3[aa]


@lru_cache()
class KideraFactor4(PeptideFeature):
    name = "Kidera_factor_4"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_4

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_4[aa]


@lru_cache()
class KideraFactor5(PeptideFeature):
    name = "Kidera_factor_5"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_5

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_5[aa]


@lru_cache()
class KideraFactor6(PeptideFeature):
    name = "Kidera_factor_6"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_6

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_6[aa]


@lru_cache()
class KideraFactor7(PeptideFeature):
    name = "Kidera_factor_7"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_7

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_7[aa]


@lru_cache()
class KideraFactor8(PeptideFeature):
    name = "Kidera_factor_8"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_8

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_8[aa]


@lru_cache()
class KideraFactor9(PeptideFeature):
    name = "Kidera_factor_9"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_9

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_9[aa]


@lru_cache()
class KideraFactor10(PeptideFeature):
    name = "Kidera_factor_10"

    @property
    @lru_cache()
    def values(self) -> dict:
        return KIDERA_FACTOR_10

    def _calculate(self, aa: str) -> float:
        return KIDERA_FACTOR_10[aa]


features_map = {
    "charge": Charge(),
    "hydrophob": Hydrophobicity(),
    "hydrophil": Hydrophilicity(),
    "isoelectric": IsoelectricPoint(),
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
    "kidera1": KideraFactor1(),
    "kidera2": KideraFactor2(),
    "kidera3": KideraFactor3(),
    "kidera4": KideraFactor4(),
    "kidera5": KideraFactor5(),
    "kidera6": KideraFactor6(),
    "kidera7": KideraFactor7(),
    "kidera8": KideraFactor8(),
    "kidera9": KideraFactor9(),
    "kidera10": KideraFactor10(),
}

operators_map = {
    "prod": ProductOperator(),
    "absdiff": AbsDifferenceOperator(),
    "diff": DifferenceOperator(),
    "layer": LayeredOperator(),
    "best": "best",
}


def parse_features(string):
    """Return a list of peptide feature objects based on a string of feature names.

    Parameters
    ----------
    string : str
        A string containing comma separated feature names.

    Returns
    -------
    list
        A list of peptide feature objects.

    Raises
    ------
    e
        Raises a value error if an unknown feature is encountered.
    """
    names = [name.lower().strip() for name in string.split(",")]
    try:
        return [features_map[name] for name in names]
    except ValueError as e:
        print("Unkown feature encountered")
        raise e


def parse_operator(string):
    """Return an operator object based on the input name.

    Parameters
    ----------
    string : str
        A string descriptor of the desired operator.

    Returns
    -------
    src.bio.operator.Operator
        The operator object matching the input name.
    """
    return operators_map[string]
