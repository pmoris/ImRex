""" Preprocess VDJdb data files in order to extract CDR3 and epitope sequence pairs.
    Turns raw VDJdb data from (../raw/vdjdb) into cleaned data ready to be analyzed (saved in ../processed).
"""

import argparse
import logging
import pandas as pd

from pathlib import Path

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
    choices=["both", "TRA", "TRB"],
    default="both",
    help='Specify which TCR chain will be extracted: "TRA", "TRB" or "both" (default)',
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


def filter_vdjdb(
    input: str, output: str, tcr_chain: str, species: str, mhc: str, hla: str
):
    """[summary]

    Parameters
    ----------
    input : str
        Filepath to the VDJdb data file, should be located in "./data/raw/vdjdb/".
    output : str
        Filepath where output should be saved, preferably "./data/processed/".
    tcr_chain : str
        Specify which TCR chain will be extracted: "TRA", "TRB" or "both" (default).
    species : str
        Specify which TCR host species will be extracted: "human" (default), "mouse", "macaque" or "all".
    mhc : str
        Specify which MHC type will be extracted: "all" (default), "MHCI" or "MHCII".
    hla : str
        Specify which HLA type will be extracted: "all" (default) or a prefix such as "HLA-A".
    """

    # initialise logger
    logger = logging.getLogger(__name__)
    logger.info("Filtering VDJdb data")

    # log the arguments that were used to call the script
    logger.info(f"Command line call: {args}")

    # read in VDJdb file (should be tab-separated)
    df = pd.read_csv(input, sep="\t")
    logger.info(f"VDJdb dataset size is {df.shape}")

    # filter rows on TCR chain
    if tcr_chain == "both":
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
        df = df.loc[df["MHC A"].str.startswith("HLA-A")]
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

    df.to_csv(output, index=False, sep=";")
    logger.info(f"Saved processed dataset to {output}.")


if __name__ == "__main__":

    input_file = Path(args.input)
    output_file = Path(args.output)
    log_file = output_file.with_suffix(".log")

    log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(filename=log_file, level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    # project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    # load_dotenv(find_dotenv())

    filter_vdjdb(
        input_file, output_file, args.tcr_chain, args.species, args.mhc, args.hla
    )
