import argparse
from collections import Counter
import logging
from pathlib import Path
import random

import pandas as pd

from src.definitions.amino_acid_properties import AMINO_ACIDS


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
        "-w",
        "--weighted",
        dest="weighted",
        action="store_true",
        help="Use weighted amino acid frequencies when generating decoys.",
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


def create_decoy_epitopes(input: str, weighted: bool = False):  # noqa: A002
    logger = logging.getLogger(__name__)

    # read in data
    df = pd.read_csv(input, sep=";")

    # get unique list of epitopes
    epitope_list = df["antigen.epitope"].unique()

    # create a dictionary mapping each epitope to a decoy one
    if weighted:
        logger.info(
            "Creating decoys using the amino acid frequencies in the original unique epitope list."
        )

        aa_distribution = dict(Counter("".join(epitope_list)))
        total = sum(aa_distribution.values())
        aa_distribution = {aa: freq / total for aa, freq in aa_distribution.items()}

        decoy_weighted = lambda epitope: "".join(
            random.choices(
                list(aa_distribution.keys()),
                list(aa_distribution.values()),
                k=len(epitope),
            )
        )
        decoy_dict = {epitope: decoy_weighted(epitope) for epitope in epitope_list}

    else:
        logger.info(
            "Creating decoys by sampling uniformly from the amino acid alphabet."
        )

        decoy = lambda epitope: "".join(random.sample(AMINO_ACIDS, len(epitope)))
        decoy_dict = {epitope: decoy(epitope) for epitope in epitope_list}

    # convert every epitope in the dataframe to its decoy
    df["antigen.epitope"] = df["antigen.epitope"].apply(
        lambda epitope: decoy_dict[epitope]
    )

    return df


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

    # log the arguments that were used to call the script
    logger.info(f"Command line call: {args}")

    # convert epitopes to decoys
    df = create_decoy_epitopes(input=input_file, weighted=args.weighted)

    # save output
    df.to_csv(output_file, index=False, sep=";")
    logger.info(f"Saved dataset with decoy epitopes to {output_file.resolve()}.")
