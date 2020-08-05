""" Preprocess immuneCODE data files in order to extract CDR3 and epitope sequence pairs.

Turns raw data from (../raw/immuneCODE) into cleaned data ready to be analyzed (saved in ../interim/immuneCODE).
"""

import argparse
import json
import logging
from pathlib import Path
import string


from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd

import pyteomics.parser
import seaborn as sns


def create_parser():
    parser = argparse.ArgumentParser(
        description="Script to extract CDR3-epitope sequence pairs from immuneCODE files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=str,
        required=True,
        help="Filepath to input immuneCODE file (peptide-detail.csv).",
    )
    parser.add_argument(
        "--length-restriction",
        dest="length_restriction",
        nargs="+",
        type=str,
        default=None,
        help="Specify the sequence length restriction. Format: cdr3-min cdr3-max epitope-min epitope-max. E.g. '10 20 8 11'. Do not use quotes around this argument.",
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


def filter_immunecode(  # noqa: C901
    input: str, length_restriction: list = None,  # noqa: A002
):
    """Filter relevant CDR3-epitope pairs from immuneCODE peptide-details file and returns a dataframe.

    Parameters
    ----------
    input : str
        Filepath to the immuneCODE data file, should be located in "./data/raw/immunecode-adaptive/ImmuneCODE-Release001.1".
    length_restriction: list
        Specify the sequence length restrictions. Format: cdr3-min cdr3-max epitope-min epitope-max. E.g. '10 20 8 11', by default None
    Returns
    -------
    DataFrame
        The filtered immuneCODE DataFrame.
    """
    # setup logger
    logger = logging.getLogger(__name__)

    logger.info(f"Filtering immuneCODE data found at {input.resolve()}")

    # read in the immuneCODE peptide-detail.csv (should be comma-separated)
    df = pd.read_csv(input, sep=",")
    logger.info(f"immuneCODE dataset size is {df.shape}")

    # assert that the immuneCODE file has the expected headers
    assert df.columns.tolist() == [
        "TCR BioIdentity",
        "TCR Nucleotide Sequence",
        "Experiment",
        "ORF Coverage",
        "Amino Acids",
        "Start Index in Genome",
        "End Index in Genome",
    ], f"Column headers of {input.resolve()} do not match the expected header names."

    # filter on rows with a single assigned peptide
    df = df[df["Amino Acids"].str.split(",").str.len() == 1]
    logger.info(
        f"Filtered down to {df.shape[0]} entries with a single assigned epitope..."
    )

    # rename column
    df = df.rename(columns={"Amino Acids": "antigen.epitope"})

    # split into CDR3, V and J gene
    df[["cdr3", "TRBV_gene", "TRBJ_gene"]] = df["TCR BioIdentity"].str.split(
        "+", expand=True
    )

    # filter on V/J genes that have a completely known allele name
    df = df[df["TRBV_gene"].str.match("TCRBV\d\d-\d\d")]
    df = df[df["TRBJ_gene"].str.match("TCRBJ.+")]
    logger.info(
        f"Filtered on full V/J gene names, resulting in {df.shape[0]} remaining entries..."
    )

    # remove "C" from V/J genes: TCRBV07-02 => TRBV07-02
    df["TRBV_gene"] = df["TRBV_gene"].str.replace("C", "")
    df["TRBJ_gene"] = df["TRBJ_gene"].str.replace("C", "")

    # filter on full cdr3 sequences (no *)
    df = df[df["cdr3"].str.match("^C\w+F$")]
    logger.info(
        f"Filtered on full CDR3 sequences, resulting in {df.shape[0]} remaining entries..."
    )

    # Filter on sequence length
    if length_restriction:
        length_restriction = [int(n) for n in length_restriction]
        assert_length_restrictions(length_restriction)
        cdr3_min, cdr3_max, epitope_min, epitope_max = length_restriction

        df = df.loc[
            (df["cdr3"].str.len() >= cdr3_min)
            & (df["cdr3"].str.len() <= cdr3_max)
            & (df["antigen.epitope"].str.len() >= epitope_min)
            & (df["antigen.epitope"].str.len() <= epitope_max)
        ]
        logger.info(
            f"Filtered on provided length restrictions ({cdr3_min}<=CDR3<={cdr3_max}, {epitope_min}<=epitope<={epitope_max}), resulting in {df.shape[0]} remaining entries..."
        )
    else:
        logger.info("Not filtering on sequence length...")

    # remove duplicates CDR3-epitope sequence pairs
    unique_columns = ["cdr3", "antigen.epitope"]
    pre_duplicate_count = df.shape[0]
    df = df.drop_duplicates(unique_columns)
    logger.info(
        f"Removing {pre_duplicate_count-df.shape[0]} duplicate CDR3-epitope pairs."
    )
    logger.info(f"Filtered dataset contains {df.shape[0]} entries.")

    return df


def visualise_data(df: pd.DataFrame, plot_file):
    sns.set_palette("Set1")
    sns.set_style("darkgrid")
    plt.rcParams["patch.edgecolor"] = "0.2"
    plt.rcParams["patch.linewidth"] = ".8"
    # plt.rc('patch', edgecolor="0.2", linewidth=.8, alpha=.75)
    # plt.rcParams['patch.alpha'] = '.75'
    # plt.rcParams['figure.figsize'] = (12, 8)

    # set font style
    plt.rcParams.update({"font.size": 11})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams[
        "font.sans-serif"
    ] = "Source Sans Pro"  # ['Fira Sans', 'Source Sans Pro']
    font = {"weight": "normal"}  # ,'size'   : 22}

    pal = sns.color_palette("BuGn_r")

    fig = plt.figure(constrained_layout=True, dpi=200, figsize=(12, 16))
    gs = GridSpec(3, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax5 = fig.add_subplot(gs[1, :])
    ax8 = fig.add_subplot(gs[2, :])

    df["cdr3.len"] = df["cdr3"].str.len()
    df["antigen.epitope.len"] = df["antigen.epitope"].str.len()

    # cdr3 length
    df["cdr3.len"] = df["cdr3"].str.len()
    sns.countplot(
        x="cdr3.len", data=df, alpha=0.8, palette=pal, ax=ax1,
    )
    # epitope length
    df["antigen.epitope.len"] = df["antigen.epitope"].str.len()
    sns.countplot(
        x="antigen.epitope.len", data=df, alpha=0.8, palette=pal, ax=ax2,
    )

    [i.set_ylabel("Number of TCR-epitope pairs") for i in [ax1, ax2]]
    # [i.set(yscale="log") for i in [ax1,ax2,ax3,ax4]]

    ax2.set_title("Epitope length distribution")
    ax2.set_xlabel("Epitope length")

    ax1.set_title("CDR3 length distribution")
    ax1.set_xlabel("CDR3 length")

    # epitope distribution
    sns.countplot(
        x="antigen.epitope",
        order=df["antigen.epitope"].value_counts().iloc[:30].index,
        data=df,
        alpha=0.8,
        palette=pal,
        ax=ax5,
    )

    n_unique_epitopes = df["antigen.epitope"].unique().shape[0]

    ax5.tick_params(labelrotation=90)
    ax5.set_title(
        f"Epitope distribution (top 30 most abundant out of {n_unique_epitopes})"
    )
    ax5.set_ylabel("Number of TCR-epitope pairs")
    ax5.set_xlabel("Epitope")

    # cdr3 counts
    df["cdr3_count"] = df.groupby("cdr3")["cdr3"].transform("size")
    # x = df.groupby("cdr3").size().sort_values(ascending=False)
    ax8 = sns.countplot(
        x="cdr3_count",
        data=df.sort_values(by="cdr3_count", ascending=False),
        alpha=0.8,
        palette=pal,
    )
    ax8.set_ylabel("Count (log scale)")
    ax8.set_xlabel("Number of epitope partners per CDR3 sequence")
    ax8.set_title("CDR3 cross-reactivity")
    # ax8.set_xlabel('Number of occurrences of a given CDR3 in the dataset')
    # ax8.set_title('CDR3 promiscuity')
    ax8.set_yscale("log")
    ax8.set_ylim((10 ** 0, 10 ** 5))

    for n, ax in enumerate(fig.axes):
        ax.text(
            -0.07,
            1.1,
            string.ascii_uppercase[n],
            transform=ax.transAxes,
            size=20,
            weight="bold",
        )

    [
        plt.setp(i.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        for i in [ax5]
    ]

    plt.savefig(plot_file)


def is_amino_acid_sequence(peptide: str):
    """Check whether a sequence contains only valid modX labels for the 20 standard amino acids.

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


def assert_length_restrictions(length_restriction_list: list):
    """Check whether number of provided length restrictions is possible.

    Parameters
    ----------
    length_restriction_list : list
        A list of length restrictions passed from the command line arguments.
    """
    assert (
        len(length_restriction_list) == 4
    ), f"Expected four numbers for length restrictions (cdr-min cdr3-max epitope-min epitope-max), but received {len(length_restriction_list)}."


if __name__ == "__main__":

    # parse cli arguments
    args = create_parser()

    # create path objects for input and output files
    input_file = Path(args.input)
    output_file = Path(args.output)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    log_file = output_file.with_suffix(".log")
    plot_file = output_file.with_suffix(".pdf")

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
    df = filter_immunecode(
        input=input_file, length_restriction=args.length_restriction,
    )

    # create visualisation of dataset
    visualise_data(df, plot_file)

    # save output
    # extract relevant columns
    columns = ["cdr3", "antigen.epitope", "TRBV_gene", "TRBJ_gene"]
    df = df.filter(items=columns)
    df.to_csv(output_file, index=False, sep=";")
    logger.info(f"Saved processed ImmuneCODE dataset to {output_file.resolve()}.")
