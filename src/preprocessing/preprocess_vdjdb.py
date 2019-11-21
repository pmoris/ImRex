""" Preprocess VDJdb data files in order to extract CDR3 and epitope sequence pairs.
    Turns raw VDJdb data from (../raw/vdjdb) into cleaned data ready to be analyzed (saved in ../processed).
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import pyteomics.parser


def create_parser():
    parser = argparse.ArgumentParser(
        description="Script to extract CDR3-epitope sequence pairs from VDJdb files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=str,
        required=True,
        help="Filepath to input VDJdb file.",
    )
    parser.add_argument(
        "--tcr-chain",
        dest="tcr_chain",
        type=str,
        choices=["all", "TRA", "TRB"],
        default="all",
        help='Specify which TCR chain will be extracted: "TRA", "TRB" or "all" (default)',
    )
    parser.add_argument(
        "--species",
        dest="species",
        type=str,
        choices=["all", "human", "mouse", "macaque"],
        default="all",
        help='Specify which TCR host species will be extracted: "human" (default), "mouse", "macaque" or "all"',
    )
    parser.add_argument(
        "--mhc",
        dest="mhc",
        type=str,
        choices=["all", "MHCI", "MHCII"],
        default="all",
        help='Specify which MHC type will be extracted: "all" (default), "MHCI" or "MHCII"',
    )
    parser.add_argument(
        "--hla",
        dest="hla",
        type=str,
        default="all",
        help='Specify which HLA type will be extracted: "all" (default) or a prefix such as "HLA-A"',
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=str,
        required=True,
        default="",
        help='Output file"',
    )
    args = parser.parse_args()
    return args


def filter_vdjdb(
    input: str, tcr_chain: str, species: str, mhc: str, hla: str,
):
    """[summary]

    Parameters
    ----------
    input : str
        Filepath to the VDJdb data file, should be located in "./data/raw/vdjdb/".
    tcr_chain : str
        Specify which TCR chain will be extracted: "TRA", "TRB" or "all" (default).
    species : str
        Specify which TCR host species will be extracted: "human" (default), "mouse", "macaque" or "all".
    mhc : str
        Specify which MHC type will be extracted: "all" (default), "MHCI" or "MHCII".
    hla : str
        Specify which HLA type will be extracted: "all" (default) or a prefix such as "HLA-A".
    """

    # setup logger
    logger = logging.getLogger(__name__)

    logger.info(f"Filtering VDJdb data found at {input.resolve()}")

    # read in VDJdb file (should be tab-separated)
    df = pd.read_csv(input, sep="\t")
    logger.info(f"VDJdb dataset size is {df.shape}")

    # assert that the VDJdb file has the expected headers
    assert_columns(df, input)

    # filter rows on TCR chain
    if tcr_chain == "all":
        logger.info("Not filtering on TCR chain...")
    else:
        df = df.loc[df["Gene"] == tcr_chain]
        logger.info(f"Filtered down to {df.shape[0]} {tcr_chain} entries...")

    # filter rows on host species
    species_dict = {
        "human": "HomoSapiens",
        "mouse": "MusMusculus",
        "macaque": "MacacaMulatta",
    }
    if species == "all":
        logger.info("Not filtering on TCR host species...")
    else:
        df = df.loc[df["Species"] == species_dict[species]]
        logger.info(f"Filtered down to {df.shape[0]} {species} entries...")

    # filter rows on MHC type
    if mhc == "all":
        logger.info("Not filtering on MHC class...")
    else:
        df = df.loc[df["MHC class"] == mhc]
        logger.info(f"Filtered down to {df.shape[0]} {mhc} entries...")

    # filter rows on HLA type
    if hla == "all":
        logger.info("Not filtering on HLA type...")
    else:
        df = df.loc[df["MHC A"].str.startswith(hla)]
        logger.info(f"Filtered down to {df.shape[0]} {hla}* entries...")

    # extract CDR3 and antigen sequence columns
    columns = ["CDR3", "Epitope"]
    df = df.filter(items=columns)

    # remove duplicates
    pre_duplicate_count = df.shape[0]
    df = df.drop_duplicates(columns)
    logger.info(
        f"Removing {pre_duplicate_count-df.shape[0]} duplicate CDR3-epitope pairs."
    )
    logger.info(f"Filtered dataset contains {df.shape[0]} entries.")

    # check if entries are valid amino acid sequences
    # is_amino_acid_sequence_vectorized = np.vectorize(is_amino_acid_sequence)
    # is_amino_acid_sequence_vectorized(df.unstack().values)
    mask = df.applymap(lambda x: is_amino_acid_sequence(x)).all(axis=1)
    if not mask.all():
        logger.warning("Removing the following invalid amino acid sequences...")
        [logger.warning(row) for row in df.loc[~mask].values]
        df = df.loc[mask]
        logger.info(f"Filtered down to {df.shape[0]} entries...")

    return df


def is_amino_acid_sequence(peptide: str):
    """Checks whether a sequence contains only valid modX labels for the 20 standard amino acids.
    std_amino_acids = ['Q', 'W', 'E', 'R', 'T', 'Y', 'I', 'P', 'A', 'S',
                   'D', 'F', 'G', 'H', 'K', 'L', 'C', 'V', 'N', 'M']

    Parameters
    ----------
    peptide : str
        A peptide sequence, either TCR or epitope.

    Returns
    -------
    bool
        True or False depending on whether the sequence contains only valid characters or not.
    """
    return all(aa in pyteomics.parser.std_amino_acids for aa in peptide)


def assert_columns(df, input):
    assert df.columns.tolist() == [
        "complex.id",
        "Gene",
        "CDR3",
        "V",
        "J",
        "Species",
        "MHC A",
        "MHC B",
        "MHC class",
        "Epitope",
        "Epitope gene",
        "Epitope species",
        "Reference",
        "Method",
        "Meta",
        "CDR3fix",
        "Score",
    ], f"Column headers of {input.resolve()} do not match the expected names."


if __name__ == "__main__":

    # parse cli arguments
    args = create_parser()

    # create path objects for input and output files
    input_file = Path(args.input)
    output_file = Path(args.output)
    log_file = output_file.with_suffix(".log")

    # file logger
    log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(filename=log_file, level=logging.INFO, format=log_fmt)
    logger = logging.getLogger(__name__)

    # console logger
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter(log_fmt))
    logger.addHandler(console)

    # not used in this stub but often useful for finding various files
    # project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    # load_dotenv(find_dotenv())

    # log the arguments that were used to call the script
    logger.info(f"Command line call: {args}")

    # preprocess vdjdb files based on passed arguments
    df = filter_vdjdb(input_file, args.tcr_chain, args.species, args.mhc, args.hla,)

    # save output
    df.to_csv(output_file, index=False, sep=";")
    logger.info(f"Saved processed dataset to {output_file.resolve()}.")
