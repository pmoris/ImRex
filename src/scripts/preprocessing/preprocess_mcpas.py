""" Preprocess immuneCODE data files in order to extract CDR3 and epitope sequence pairs.

Turns raw data from (../raw/immuneCODE) into cleaned data ready to be analyzed (saved in ../interim/immuneCODE).
"""

import argparse
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
        description="Script to extract CDR3 (beta)-epitope sequence pairs from McPAS.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=str,
        required=True,
        help="Filepath to input McPAS file.",
    )
    parser.add_argument(
        "--vdjdb_dataset",
        dest="vdjdb_dataset",
        type=str,
        required=True,
        help="Filepath to (preprocessed) VDJdb data used for comparing shared and unique epitopes.",
    )
    parser.add_argument(
        "--species",
        dest="species",
        type=str,
        choices=["all", "human", "mouse", "macaque"],
        default="human",
        help='Specify which TCR host species will be extracted: "human" (default), "mouse", "macaque" or "all"',
    )
    parser.add_argument(
        "--mhc",
        dest="mhc",
        type=str,
        choices=["all", "MHCI", "MHCII"],
        default="MHCI",
        help='Specify which MHC type will be extracted: "all", "MHCI" (default) or "MHCII"',
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


def filter_mcpas(  # noqa: C901
    input: str,
    species: str = "human",
    mhc: str = "MHCI",
    length_restriction: list = None,  # noqa: A002
):
    """Filter relevant CDR3-epitope pairs from McPAS peptide-details file and returns a dataframe.

    Parameters
    ----------
    input : str
        Filepath to the McPAS data file, should be located in "./data/raw/McPAS-TCR.csv".
    species : str
        Specify which TCR host species will be extracted: "human" (default), "mouse", "macaque" or "all".
    mhc : str
        Specify which MHC type will be extracted: "all" (default), "MHCI" or "MHCII".
    length_restriction: list
        Specify the sequence length restrictions. Format: cdr3-min cdr3-max epitope-min epitope-max. E.g. '10 20 8 11', by default None
    Returns
    -------
    DataFrame
        The filtered McPAS DataFrame.
    """
    # setup logger
    logger = logging.getLogger(__name__)

    logger.info(f"Filtering McPAS data found at {input.resolve()}")

    # read in the McPAS peptide-detail.csv (should be comma-separated)
    df_mcpas = pd.read_csv(input, sep=",", encoding="latin")
    logger.info(f"McPAS dataset size is {df_mcpas.shape}")

    # assert that the McPAS file has the expected headers
    print(df_mcpas.columns.tolist())
    assert df_mcpas.columns.tolist() == [
        "CDR3.alpha.aa",
        "CDR3.beta.aa",
        "Species",
        "Category",
        "Pathology",
        "Pathology.Mesh.ID",
        "Additional.study.details",
        "Antigen.identification.method",
        "Single.cell",
        "NGS",
        "Antigen.protein",
        "Protein.ID",
        "Epitope.peptide",
        "Epitope.ID",
        "MHC",
        "Tissue",
        "T.Cell.Type",
        "T.cell.characteristics",
        "CDR3.alpha.nt",
        "TRAV",
        "TRAJ",
        "TRBV",
        "TRBD",
        "TRBJ",
        "Reconstructed.J.annotation",
        "CDR3.beta.nt",
        "Mouse.strain",
        "PubMed.ID",
        "Remarks",
        "CDR3.beta.clean",
    ], f"Column headers of {input.resolve()} do not match the expected header names."

    # rename columns
    df_mcpas = df_mcpas.rename(
        columns={"CDR3.beta.clean": "cdr3", "Epitope.peptide": "antigen.epitope"}
    )

    # create merged column
    df_mcpas["merged"] = df_mcpas["cdr3"] + "-" + df_mcpas["antigen.epitope"]

    # filter rows on host species
    species_dict = {
        "human": "Human",
        "mouse": "Mouse",
    }
    if species == "all":
        logger.info("Not filtering on TCR host species...")
    else:
        df_mcpas = df_mcpas.loc[df_mcpas["Species"] == species_dict[species]]
        logger.info(f"Filtered down to {df_mcpas.shape[0]} {species} entries...")

    # filter rows on MHC type
    mhc_dict = {
        "MHCI": "CD8",
        "MHCII": "CD4",
    }
    if mhc == "all":
        logger.info("Not filtering on MHC class / T-cell type...")
    else:
        df_mcpas = df_mcpas.loc[df_mcpas["T.Cell.Type"] == mhc_dict[mhc]]
        logger.info(f"Filtered down to {df_mcpas.shape[0]} {mhc} entries...")

    # filter on tetramers binding and/or T cell stimulation to avoid bystander effects
    df_mcpas[df_mcpas["Antigen.identification.method"].isin([1.0, 2.1])]
    logger.info(
        f"Filtered down to {df_mcpas.shape[0]} entries determined by tetramer assay or stimulation with a peptide..."
    )

    # remove entries with remarks
    df_mcpas = df_mcpas[pd.isnull(df_mcpas["Remarks"])]
    logger.info(f"Filtered down to {df_mcpas.shape[0]} entries without remarks...")

    # Filter on sequence length
    if length_restriction:
        length_restriction = [int(n) for n in length_restriction]
        assert_length_restrictions(length_restriction)
        cdr3_min, cdr3_max, epitope_min, epitope_max = length_restriction

        df_mcpas = df_mcpas.loc[
            (df_mcpas["cdr3"].str.len() >= cdr3_min)
            & (df_mcpas["cdr3"].str.len() <= cdr3_max)
            & (df_mcpas["antigen.epitope"].str.len() >= epitope_min)
            & (df_mcpas["antigen.epitope"].str.len() <= epitope_max)
        ]
        logger.info(
            f"Filtered on provided length restrictions ({cdr3_min}<=CDR3<={cdr3_max}, {epitope_min}<=epitope<={epitope_max}), resulting in {df_mcpas.shape[0]} remaining entries..."
        )
    else:
        logger.info("Not filtering on sequence length...")

    # remove duplicates CDR3-epitope sequence pairs
    unique_columns = ["cdr3", "antigen.epitope"]
    pre_duplicate_count = df_mcpas.shape[0]
    df_mcpas = df_mcpas.drop_duplicates(unique_columns)
    logger.info(
        f"Removing {pre_duplicate_count-df_mcpas.shape[0]} duplicate CDR3-epitope pairs."
    )

    # remove missing data
    pre_missing_count = df_mcpas.shape[0]
    df_mcpas = df_mcpas.dropna(subset=unique_columns)
    logger.info(
        f"Removing {pre_missing_count-df_mcpas.shape[0]} CDR3-epitope pairs with missing values."
    )

    logger.info(f"Filtered dataset contains {df_mcpas.shape[0]} entries.")

    return df_mcpas


def filter_unique_epitopes(df_mcpas: pd.DataFrame, df_vdjdb: pd.DataFrame):
    """Filter on epitopes in the McPAS dataset that do not occur in the VDJdb dataset.

    Parameters
    ----------
    df_mcpas : pd.DataFrame
        Preprocessed McPAS dataframe.
    vdjdb_dataset: pd.DataFrame
        Preprocessed VDJdb dataframe.
    Returns
    -------
    DataFrame
        The filtered McPAS DataFrame.
    """

    logger.info(
        f"Extracting epitopes that are NOT present in the supplied VDJdb dataset {Path(args.vdjdb_dataset)} entries."
    )

    unique = set(df_mcpas["antigen.epitope"].unique()) - set(
        df_vdjdb["antigen.epitope"].unique()
    )
    df_mcpas_unique = df_mcpas[df_mcpas["antigen.epitope"].isin(unique)]
    logger.info(
        f"Filtered unique epitopes dataset contains {df_mcpas_unique.shape[0]} entries."
    )

    return df_mcpas_unique


def filter_shared_epitopes(df_mcpas: pd.DataFrame, df_vdjdb: pd.DataFrame):
    """Filter on epitopes in the McPAS dataset that do also occur in the VDJdb dataset.

    Parameters
    ----------
    df_mcpas : pd.DataFrame
        Preprocessed McPAS dataframe.
    vdjdb_dataset: pd.DataFrame
        Preprocessed VDJdb dataframe.
    Returns
    -------
    DataFrame
        The filtered McPAS DataFrame.
    """

    logger.info(
        f"Extracting epitopes that ARE present in the supplied VDJdb dataset {Path(args.vdjdb_dataset)} entries."
    )

    shared = set(df_mcpas["antigen.epitope"].unique()).intersection(
        set(df_vdjdb["antigen.epitope"].unique())
    )
    df_mcpas_shared = df_mcpas[df_mcpas["antigen.epitope"].isin(shared)]

    logger.info(
        f"McPAS shares {len(shared)} epitopes with VDJdb dataset, resulting in {df_mcpas_shared.shape[0]} entries."
    )

    # remove identical entries to VDJdb
    df_mcpas_shared = (
        df_mcpas_shared[["cdr3", "antigen.epitope"]]
        .merge(df_vdjdb[["cdr3", "antigen.epitope"]], indicator=True, how="left")
        .loc[lambda x: x["_merge"] == "left_only"]
    )
    logger.info(
        f"Removing CDR3-epitope pairs that are identical to those in the supplied VDJdb dataset {Path(args.vdjdb_dataset)} entries."
    )

    logger.info(
        f"Filtered shared epitopes dataset contains {df_mcpas_shared.shape[0]} unique entries."
    )

    return df_mcpas_shared


def visualise_data(df_mcpas: pd.DataFrame, plot_file):
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

    df_mcpas["cdr3.len"] = df_mcpas["cdr3"].str.len()
    df_mcpas["antigen.epitope.len"] = df_mcpas["antigen.epitope"].str.len()

    # cdr3 length
    df_mcpas["cdr3.len"] = df_mcpas["cdr3"].str.len()
    sns.countplot(
        x="cdr3.len", data=df_mcpas, alpha=0.8, palette=pal, ax=ax1,
    )
    # epitope length
    df_mcpas["antigen.epitope.len"] = df_mcpas["antigen.epitope"].str.len()
    sns.countplot(
        x="antigen.epitope.len", data=df_mcpas, alpha=0.8, palette=pal, ax=ax2,
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
        order=df_mcpas["antigen.epitope"].value_counts().iloc[:30].index,
        data=df_mcpas,
        alpha=0.8,
        palette=pal,
        ax=ax5,
    )

    n_unique_epitopes = df_mcpas["antigen.epitope"].unique().shape[0]

    ax5.tick_params(labelrotation=90)
    ax5.set_title(
        f"Epitope distribution (top 30 most abundant out of {n_unique_epitopes})"
    )
    ax5.set_ylabel("Number of TCR-epitope pairs")
    ax5.set_xlabel("Epitope")

    # cdr3 counts
    df_mcpas["cdr3_count"] = df_mcpas.groupby("cdr3")["cdr3"].transform("size")
    # x = df_mcpas.groupby("cdr3").size().sort_values(ascending=False)
    ax8 = sns.countplot(
        x="cdr3_count",
        data=df_mcpas.sort_values(by="cdr3_count", ascending=False),
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
    vdjdb_dataset = Path(args.vdjdb_dataset)
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
    df_mcpas = filter_mcpas(
        input=input_file,
        species=args.species,
        mhc=args.mhc,
        length_restriction=args.length_restriction,
    )

    # create visualisation of dataset
    visualise_data(df_mcpas, plot_file)

    # save output
    # extract relevant columns
    columns = ["cdr3", "antigen.epitope", "TRBV_gene", "TRBJ_gene"]
    df_mcpas = df_mcpas.filter(items=columns)
    df_mcpas.to_csv(output_file, index=False, sep=";")
    logger.info(f"Saved processed McPAS dataset to {output_file.resolve()}.")

    # read in the preprocessed VDJdb dataset used for comparing epitopes
    df_vdjdb = pd.read_csv(vdjdb_dataset, sep=";")
    logger.info(f"VDJdb dataset size is {df_vdjdb.shape}")

    # repeat for unique epitopes
    df_mcpas_unique = filter_unique_epitopes(df_mcpas.copy(), df_vdjdb)

    unique_output_file = output_file.parent / (output_file.stem + "-unique.csv")
    unique_plot_file = unique_output_file.with_suffix(".pdf")

    visualise_data(df_mcpas_unique, unique_plot_file)

    df_mcpas_unique.to_csv(unique_output_file, index=False, sep=";")

    logger.info(
        f"Saved processed unique epitope McPAS dataset to {unique_output_file.resolve()}."
    )

    # repeat for shared epitopes
    df_mcpas_shared = filter_shared_epitopes(df_mcpas.copy(), df_vdjdb)

    shared_output_file = output_file.parent / (output_file.stem + "-shared.csv")
    shared_plot_file = shared_output_file.with_suffix(".pdf")

    visualise_data(df_mcpas_shared, shared_plot_file)

    df_mcpas_shared.to_csv(shared_output_file, index=False, sep=";")

    logger.info(
        f"Saved processed shared epitope McPAS dataset to {shared_output_file.resolve()}."
    )
