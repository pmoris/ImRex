import numpy as np
import pandas as pd

# list of valid amino acid characters
# taken from biopython Bio.Alphabet.IUPAC.IUPACProtein.letters
# https://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.IUPACProtein-class.html
# which might be depracated at some point in the future.

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# PH value, required for certain physiochemical properties
PH = 7

# Kidera factors
# https://github.com/vadimnazarov/kidera-atchley
# http://www.jbsdonline.com/mc_images/category/4308/14-scheraga_jbsd_28_4.pdf
_kidera_pd = pd.DataFrame.from_records(
    list(
        map(
            lambda x: list(map(float, x.split(","))),
            "-1.56,-1.67,-0.97,-0.27,-0.93,-0.78,-0.2,-0.08,0.21,-0.48;0.22,1.27,1.37,1.87,-1.7,0.46,0.92,-0.39,0.23,0.93;1.14,-0.07,-0.12,0.81,0.18,0.37,-0.09,1.23,1.1,-1.73;0.58,-0.22,-1.58,0.81,-0.92,0.15,-1.52,0.47,0.76,0.7;0.12,-0.89,0.45,-1.05,-0.71,2.41,1.52,-0.69,1.13,1.1;-0.47,0.24,0.07,1.1,1.1,0.59,0.84,-0.71,-0.03,-2.33;-1.45,0.19,-1.61,1.17,-1.31,0.4,0.04,0.38,-0.35,-0.12;1.46,-1.96,-0.23,-0.16,0.1,-0.11,1.32,2.36,-1.66,0.46;-0.41,0.52,-0.28,0.28,1.61,1.01,-1.85,0.47,1.13,1.63;-0.73,-0.16,1.79,-0.77,-0.54,0.03,-0.83,0.51,0.66,-1.78;-1.04,0,-0.24,-1.1,-0.55,-2.05,0.96,-0.76,0.45,0.93;-0.34,0.82,-0.23,1.7,1.54,-1.62,1.15,-0.08,-0.48,0.6;-1.4,0.18,-0.42,-0.73,2,1.52,0.26,0.11,-1.27,0.27;-0.21,0.98,-0.36,-1.43,0.22,-0.81,0.67,1.1,1.71,-0.44;2.06,-0.33,-1.15,-0.75,0.88,-0.45,0.3,-2.3,0.74,-0.28;0.81,-1.08,0.16,0.42,-0.21,-0.43,-1.89,-1.15,-0.97,-0.23;0.26,-0.7,1.21,0.63,-0.1,0.21,0.24,-1.15,-0.56,0.19;0.3,2.1,-0.72,-1.57,-1.16,0.57,-0.48,-0.4,-2.3,-0.6;1.38,1.48,0.8,-0.56,0,-0.68,-0.31,1.03,-0.05,0.53;-0.74,-0.71,2.04,-0.4,0.5,-0.81,-1.07,0.06,-0.46,0.65".split(
                ";"
            ),
        )
    ),
    index=[
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ],
    columns=list(map(lambda x: "f" + str(x), range(1, 11))),
)

KIDERA_FACTOR_1 = _kidera_pd.loc[:, "f1"].to_dict()
KIDERA_FACTOR_2 = _kidera_pd.loc[:, "f2"].to_dict()
KIDERA_FACTOR_3 = _kidera_pd.loc[:, "f3"].to_dict()
KIDERA_FACTOR_4 = _kidera_pd.loc[:, "f4"].to_dict()
KIDERA_FACTOR_5 = _kidera_pd.loc[:, "f5"].to_dict()
KIDERA_FACTOR_6 = _kidera_pd.loc[:, "f6"].to_dict()
KIDERA_FACTOR_7 = _kidera_pd.loc[:, "f7"].to_dict()
KIDERA_FACTOR_8 = _kidera_pd.loc[:, "f8"].to_dict()
KIDERA_FACTOR_9 = _kidera_pd.loc[:, "f9"].to_dict()
KIDERA_FACTOR_10 = _kidera_pd.loc[:, "f10"].to_dict()

# Atchley factors
# Described in Atchley et al. 2005 - https://www.pnas.org/content/102/18/6395.full (https://doi.org/10.1073/pnas.0408677102)
# Values taken from the R package HDMD: https://github.com/cran/HDMD/blob/master/R/HDMD_package.R
# Factor Scores as computed by SAS Factor Analysis of 54 Amino Acid Indices
# Rows are alphabetized Amino Acids and the 5 columns are factors where
# Factor1 (PAH): Polarity, Accessibility, Hydrophobicity
# Factor2 (PSS): Propensity for Secondary Structure
# Factor3 (MS) : Molecular Size
# Factor4 (CC): Codon Composition
# Factor5 (EC): Electrostatic Charge.

_atchley_array = np.array(
    [
        -0.59145974,
        -1.30209266,
        -0.7330651,
        1.5703918,
        -0.14550842,
        -1.34267179,
        0.46542300,
        -0.8620345,
        -1.0200786,
        -0.25516894,
        1.05015062,
        0.30242411,
        -3.6559147,
        -0.2590236,
        -3.24176791,
        1.35733226,
        -1.45275578,
        1.4766610,
        0.1129444,
        -0.83715681,
        -1.00610084,
        -0.59046634,
        1.8909687,
        -0.3966186,
        0.41194139,
        -0.38387987,
        1.65201497,
        1.3301017,
        1.0449765,
        2.06385566,
        0.33616543,
        -0.41662780,
        -1.6733690,
        -1.4738898,
        -0.07772917,
        -1.23936304,
        -0.54652238,
        2.1314349,
        0.3931618,
        0.81630366,
        1.83146558,
        -0.56109831,
        0.5332237,
        -0.2771101,
        1.64762794,
        -1.01895162,
        -0.98693471,
        -1.5046185,
        1.2658296,
        -0.91181195,
        -0.66312569,
        -1.52353917,
        2.2194787,
        -1.0047207,
        1.21181214,
        0.94535614,
        0.82846219,
        1.2991286,
        -0.1688162,
        0.93339498,
        0.18862522,
        2.08084151,
        -1.6283286,
        0.4207004,
        -1.39177378,
        0.93056541,
        -0.17926549,
        -3.0048731,
        -0.5025910,
        -1.85303476,
        1.53754853,
        -0.05472897,
        1.5021086,
        0.4403185,
        2.89744417,
        -0.22788299,
        1.39869991,
        -4.7596375,
        0.6701745,
        -2.64747356,
        -0.03181782,
        0.32571153,
        2.2134612,
        0.9078985,
        1.31337035,
        -1.33661279,
        -0.27854634,
        -0.5440132,
        1.2419935,
        -1.26225362,
        -0.59533918,
        0.00907760,
        0.6719274,
        -2.1275244,
        -0.18358096,
        0.25999617,
        0.82992312,
        3.0973596,
        -0.8380164,
        1.51150958,
    ]
).reshape(20, 5)

ATCHLEY_FACTOR_1 = {}
ATCHLEY_FACTOR_2 = {}
ATCHLEY_FACTOR_3 = {}
ATCHLEY_FACTOR_4 = {}
ATCHLEY_FACTOR_5 = {}

for index, aa in enumerate(AMINO_ACIDS):
    ATCHLEY_FACTOR_1[aa] = _atchley_array[index, 0]
    ATCHLEY_FACTOR_2[aa] = _atchley_array[index, 1]
    ATCHLEY_FACTOR_3[aa] = _atchley_array[index, 2]
    ATCHLEY_FACTOR_4[aa] = _atchley_array[index, 3]
    ATCHLEY_FACTOR_5[aa] = _atchley_array[index, 4]

# TCRex amino acid features

# Source: https://github.com/bittremieux/TCR-Classifier/blob/master/tcr_classifier_v2.ipynb
TCREX_BASICITY = {
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

TCREX_HYDROPHOBICITY = {
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

TCREX_HELICITY = {
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

TCREX_MUTATION_STABILITY = {
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
