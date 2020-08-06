""" Preprocess VDJdb data files in order to extract CDR3 and epitope sequence pairs.

Turns raw VDJdb data from (../raw/vdjdb) into cleaned data ready to be analyzed (saved in ../processed).
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
        description="Script to extract CDR3-epitope sequence pairs from VDJdb files and apply various filters and preprocessing steps.",
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
        help='Specify which TCR host species will be extracted: "human", "mouse", "macaque" or "all" (default)',
    )
    spurious_parser = parser.add_mutually_exclusive_group(required=False)
    spurious_parser.add_argument(
        "--drop-spurious",
        dest="drop_spurious",
        action="store_true",
        help="Indicates that spurious sequence pairs (as defined by cdr3fix: good = false) should be dropped (default behaviour is to drop them).",
    )
    spurious_parser.add_argument(
        "--keep-spurious",
        dest="drop_spurious",
        action="store_false",
        help="Indicates that spurious sequence pairs (as defined by cdr3fix: good = false) should be kept.",
    )
    parser.set_defaults(drop_spurious="True")
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
        "--remove-specific-epitope-reference",
        dest="specific_removal_epitope_reference",
        nargs="+",
        type=str,
        default=None,
        help="Specify a specific epitope and reference to remove. E.g. 'KLGGALQAK 10xgenomics'",
    )
    parser.add_argument(
        "--remove-specific-references",
        dest="specific_removal_references",
        nargs="+",
        type=str,
        default=None,
        help="Specify specific references to remove. E.g. '10xgenomics PMID#####'",
    )
    parser.add_argument(
        "--keep-specific-references",
        dest="keep_specific_references",
        nargs="+",
        type=str,
        default=None,
        help="Specify specific references to keep. E.g. '10xgenomics PMID#####'",
    )
    parser.add_argument(
        "--keep-specific-epitopes",
        dest="keep_specific_epitopes",
        nargs="+",
        type=str,
        default=None,
        help="Specify specific epitopes to keep. E.g. 'NLVPMVATV KLGGALQAK'",
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
        "--downsample",
        dest="downsample",
        nargs="+",
        type=str,
        default=None,
        help="Specify which epitopes should be downsampled. Format: epitope-seq fraction-to-drop. E.g. 'NLVPMVATV 1000 GILGFVFTL 1000'. Do not use quotes around this argument.",
    )
    parser.add_argument(
        "--terminal_only",
        dest="terminal_only",
        action="store_true",
        default=False,
        help="Only keep the two terminal amino acids on each side of the CDR3 (i.e. 4 aa in total).",
    )
    parser.add_argument(
        "--middle_only",
        dest="middle_only",
        action="store_true",
        default=False,
        help="Trim the two terminal amino acids on each side of the CDR3 (i.e. 4 aa in total).",
    )
    parser.add_argument(
        "--terminal_replaced",
        dest="terminal_replaced",
        type=str,
        choices=["X", "G"],
        default=None,
        help="Replace the two terminal amino acids on each side of the CDR3 (i.e. 4 aa in total) by an X or glycine.",
    )
    parser.add_argument(
        "--middle_replaced",
        type=str,
        choices=["X", "G"],
        default=None,
        help="Replace everything but the two terminal amino acids on each side of the CDR3 (i.e. 4 aa in total) by an X or glycine.",
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


def filter_vdjdb(  # noqa: C901
    input: str,  # noqa: A002
    tcr_chain: str = "all",
    species: str = "all",
    drop_spurious: bool = True,
    mhc: str = "all",
    hla: str = "all",
    specific_removal_epitope_reference: list = None,
    specific_removal_references: list = None,
    keep_specific_references: list = None,
    keep_specific_epitopes: list = None,
    length_restriction: list = None,
    downsample: list = None,
    terminal_only=False,
    middle_only=False,
    terminal_replaced=None,
    middle_replaced=None,
):
    """Filter relevant CDR3-epitope pairs from VDJdb files and returns a dataframe.

    Parameters
    ----------
    input : str
        Filepath to the VDJdb data file, should be located in "./data/raw/vdjdb/".
    tcr_chain : str
        Specify which TCR chain will be extracted: "TRA", "TRB" or "all" (default).
    species : str
        Specify which TCR host species will be extracted: "human" (default), "mouse", "macaque" or "all".
    drop_spurious : boolan
        Indicates whether or not the spurious sequences (as defined by cdr3fix: good = false) should be dropped (default: True, i.e. drop them).
    mhc : str
        Specify which MHC type will be extracted: "all" (default), "MHCI" or "MHCII".
    hla : str
        Specify which HLA type will be extracted: "all" (default) or a prefix such as "HLA-A".
    specific_removal_epitope_reference : list
        Specify a specific epitope and reference to remove. E.g. 'KLGGALQAK 10xgenomics'", by default None
    specific_removal_references : list
        Specify specific references to remove. E.g. '10xgenomics PMID#####', by default None
    keep_specific_references : list
        Specify specific references to keep. E.g. '10xgenomics PMID#####', by default None
    keep_specific_epitopes : list
        Specify specific epitopes to keep. E.g. 'NLVPMVATV KLGGALQAK', by default None
    length_restriction: list
        Specify the sequence length restrictions. Format: cdr3-min cdr3-max epitope-min epitope-max. E.g. '10 20 8 13', by default None
    Downsample: list
        Specify which epitopes should be downsampled. Format: epitope-seq fraction-to-drop. E.g. 'NLVPMVATV 0.84 GILGFVFTL 0.80', by default None
    terminal_only: bool
        Only keep the first and the last two amino acids of the CDR3 sequence.
    middle_only: bool
        Trim the first and the last two amino acids of the CDR3 sequence.
    terminal_replaced: str
        Replace the first and the last two amino acids of the CDR3 sequence by either an X or a G.
    middle_replaced: str,
        Replace everything but the first and the last two amino acids of the CDR3 sequence by either an X or a G.

    Returns
    -------
    DataFrame
        The filtered vdjdb DataFrame.
    """
    # setup logger
    logger = logging.getLogger(__name__)

    logger.info(f"Filtering VDJdb data found at {input.resolve()}")

    # read in VDJdb file (should be tab-separated)
    df = pd.read_csv(input, sep="\t")
    logger.info(f"VDJdb dataset size is {df.shape}")

    # assert that the VDJdb file has the expected headers
    assert_columns(df, input)

    # parse the dict/json-like methods, meta and cdr3fix columns
    df = (
        df.join(
            pd.json_normalize(
                df["method"].apply(lambda x: json.loads(r"{}".format(x)))
            ).add_prefix("method.")
        )
        .join(
            pd.json_normalize(
                df["meta"].apply(lambda x: json.loads(r"{}".format(x)))
            ).add_prefix("meta.")
        )
        .join(
            pd.json_normalize(
                df["cdr3fix"].apply(lambda x: json.loads(r"{}".format(x)))
            ).add_prefix("cdr3fix.")
        )
        .drop(["method", "meta", "cdr3fix"], axis=1)
    )

    # filter rows on TCR chain
    if tcr_chain == "all":
        logger.info("Not filtering on TCR chain...")
    else:
        # df = df.loc[df["Gene"] == tcr_chain]
        df = df.loc[df["gene"] == tcr_chain]
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
        # df = df.loc[df["Species"] == species_dict[species]]
        df = df.loc[df["species"] == species_dict[species]]
        logger.info(f"Filtered down to {df.shape[0]} {species} entries...")

    # remove spurious entries
    if drop_spurious:
        df = df.loc[df["cdr3fix.good"]]
        logger.info(f"Filtered down to {df.shape[0]} non-spurious entries...")
    else:
        logger.info("Not removing spurious entries...")

    # filter rows on MHC type
    if mhc == "all":
        logger.info("Not filtering on MHC class...")
    else:
        # df = df.loc[df["MHC class"] == mhc]
        df = df.loc[df["mhc.class"] == mhc]
        logger.info(f"Filtered down to {df.shape[0]} {mhc} entries...")

    # filter rows on HLA type
    if hla == "all":
        logger.info("Not filtering on HLA type...")
    else:
        # df = df.loc[df["MHC A"].str.startswith(hla)]
        df = df.loc[df["mhc.a"].str.startswith(hla)]
        logger.info(f"Filtered down to {df.shape[0]} {hla}* entries...")

    # Remove specifically provided epitope-reference combination
    if specific_removal_epitope_reference:
        df = df.loc[
            ~(df["antigen.epitope"] == specific_removal_epitope_reference[0])
            | ~(df["reference.id"].str.contains(specific_removal_epitope_reference[1]))
        ]
        logger.info(
            f"Removed specific entries with epitope sequence {specific_removal_epitope_reference[0]} and reference {specific_removal_epitope_reference[1]}, resulting in {df.shape[0]} remaining entries..."
        )
    else:
        logger.info(
            "Not removing any specific epitope-reference combination entries..."
        )

    # Remove specified references
    if specific_removal_references:
        df = df.loc[
            ~(df["reference.id"].str.contains("|".join(specific_removal_references)))
        ]
        # df["reference.id"].isin([specific_removal_references]) could be used if exact matching is wanted instead of containing the string
        logger.info(
            f"Removed specific entries with reference(s) {specific_removal_references}, resulting in {df.shape[0]} remaining entries..."
        )
    else:
        logger.info("Not removing any specific reference entries...")

    # Keep specified references
    if keep_specific_references:
        df = df.loc[
            (df["reference.id"].str.contains("|".join(keep_specific_references)))
        ]
        logger.info(
            f"Filtered on specific entries with reference(s) {keep_specific_references}, resulting in {df.shape[0]} remaining entries..."
        )
    else:
        logger.info("Not filtering on specific reference entries...")

    if keep_specific_epitopes:
        df = df.loc[(df["antigen.epitope"].isin(keep_specific_epitopes))]
        logger.info(
            f"Filtered on specific epitopes {keep_specific_references}, resulting in {df.shape[0]} remaining entries..."
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

    # Downsample epitopes - should happen after de-duplication!
    if downsample:
        for epitope, n in zip(downsample[::2], downsample[1::2]):
            # get indices for n randomly sampled rows
            sampled_indices = (
                df[df["antigen.epitope"] == epitope]
                .sample(n=int(n), random_state=42)
                .index
            )
            # select inverse indices to drop from dataframe
            drop_indices = df[df["antigen.epitope"] == epitope].index.difference(
                sampled_indices
            )
            # drop those from the dataframe
            df = df.drop(drop_indices)
            # check remaining rows for this epitope
            m = len(drop_indices)
            logger.info(
                f"Filtered down to {n} {epitope} observations, after removing {m}."
            )
        logger.info(f"After downsampling, there are {df.shape[0]} remaining entries...")

    # check if entries are valid amino acid sequences
    # is_amino_acid_sequence_vectorized = np.vectorize(is_amino_acid_sequence)
    # is_amino_acid_sequence_vectorized(df.unstack().values)
    mask = df[unique_columns].applymap(lambda x: is_amino_acid_sequence(x)).all(axis=1)
    if not mask.all():
        logger.warning("Removing the following invalid amino acid sequences...")
        [logger.warning(row) for row in df.loc[~mask].values]
        df = df.loc[mask]
        logger.info(f"Filtered down to {df.shape[0]} entries...")

    # replace terminal ends / middle region
    # needs to happen at the very end of the filtering steps and requires an additional deduplication round
    if terminal_only:
        pattern = r"(?P<terminalstart>^.{2})(?P<middle>.*)(?P<terminalend>.{2}$)"
        df["cdr3"] = df["cdr3"].str.replace(
            pattern, lambda m: m.group("terminalstart") + m.group("terminalend")
        )

    elif middle_only:
        pattern = r"(?P<terminalstart>^.{2})(?P<middle>.*)(?P<terminalend>.{2}$)"
        df["cdr3"] = df["cdr3"].str.replace(pattern, lambda m: m.group("middle"),)

    elif terminal_replaced:
        df["cdr3"] = df["cdr3"].str.replace("^..", terminal_replaced * 2)
        df["cdr3"] = df["cdr3"].str.replace("..$", terminal_replaced * 2)

    elif middle_replaced:
        pattern = r"(?P<terminalstart>^.{2})(?P<middle>.*)(?P<terminalend>.{2}$)"
        df["cdr3"] = df["cdr3"].str.replace(
            pattern,
            lambda m: m.group("terminalstart")
            + middle_replaced * len(m.group("middle"))
            + m.group("terminalend"),
        )

    # remove duplicates created by replacing amino acids
    if any([terminal_only, middle_only, terminal_replaced, middle_replaced,]):
        pre_duplicate_count = df.shape[0]
        df = df.drop_duplicates(unique_columns)
        logger.info(
            f"Removing {pre_duplicate_count-df.shape[0]} duplicate CDR3-epitope pairs that were introduced after trimming or replacing middle/terminal amino acids in the CDR3 sequence."
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
    gs = GridSpec(6, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, :])
    ax6 = fig.add_subplot(gs[3, :])
    ax7 = fig.add_subplot(gs[4, :])
    ax8 = fig.add_subplot(gs[5, :])

    df["cdr3.len"] = df["cdr3"].str.len()
    df["antigen.epitope.len"] = df["antigen.epitope"].str.len()

    # cdr3 length - TRB

    if not df.loc[df["mhc.class"] == "MHCI"].empty:
        sns.countplot(
            x="cdr3.len",
            hue="gene",
            data=df.loc[df["mhc.class"] == "MHCI"].sort_values(by="gene"),
            alpha=0.8,
            palette=pal,
            ax=ax1,
        )
    # epitope length - TRB
    if not df.loc[df["mhc.class"] == "MHCI"].empty:
        sns.countplot(
            x="antigen.epitope.len",
            hue="gene",
            data=df.loc[df["mhc.class"] == "MHCI"].sort_values(by="gene"),
            alpha=0.8,
            palette=pal,
            ax=ax2,
        )
    # cdr3 length - TRA
    if not df.loc[df["mhc.class"] == "MHCII"].empty:
        sns.countplot(
            x="cdr3.len",
            hue="gene",
            data=df.loc[df["mhc.class"] == "MHCII"].sort_values(by="gene"),
            alpha=0.8,
            palette=pal,
            ax=ax3,
        )
    # epitope length - TRA
    if not df.loc[df["mhc.class"] == "MHCII"].empty:
        sns.countplot(
            x="antigen.epitope.len",
            hue="gene",
            data=df.loc[df["mhc.class"] == "MHCII"].sort_values(by="gene"),
            alpha=0.8,
            palette=pal,
            ax=ax4,
        )

    [i.set_ylabel("Number of TCR-epitope pairs") for i in [ax1, ax2, ax3, ax4]]
    # [i.set(yscale="log") for i in [ax1,ax2,ax3,ax4]]
    [i.legend(title="TCR chain", loc="upper left") for i in [ax1, ax2, ax3, ax4]]

    ax2.set_title("Epitope length distribution for MHCI")
    ax4.set_title("Epitope length distribution for MHCII")
    ax2.set_xlabel("Epitope length")
    ax4.set_xlabel("Epitope length")

    ax1.set_title("CDR3 length distribution for MHCI")
    ax3.set_title("CDR3 length distribution for MHCII")
    ax1.set_xlabel("CDR3 length")
    ax3.set_xlabel("CDR3 length")

    # epitope distribution
    sns.countplot(
        x="antigen.epitope",
        hue="gene",
        order=df["antigen.epitope"].value_counts().iloc[:30].index,
        data=df.sort_values(by="gene"),
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
    ax5.legend(title="TCR chain", loc="upper right")

    # MHC distribution A
    sns.countplot(
        x="mhc.a",
        hue="gene",
        order=df["mhc.a"].value_counts().index,
        data=df.sort_values(by="gene"),
        palette=pal,
        alpha=0.8,
        ax=ax6,
    )

    ax6.tick_params(labelrotation=90)
    ax6.set_ylabel("Number of TCR-epitope pairs")
    ax6.set_xlabel("MHC alpha subunit")
    ax6.set_title("MHC alpha subunit allele distribution")
    ax6.legend(title="TCR chain", loc="upper right")

    # MHC distribution B
    sns.countplot(
        x="mhc.b",
        hue="gene",
        order=df["mhc.b"].value_counts().index,
        data=df.sort_values(by="gene"),
        palette=pal,
        alpha=0.8,
        ax=ax7,
    )

    ax7.tick_params(labelrotation=90)
    ax7.set_ylabel("Number of TCR-epitope pairs (log scale)")
    ax7.set_xlabel("MHC beta subunit")
    ax7.set_title("MHC beta subunit allele distribution")
    ax7.legend(title="TCR chain", loc="upper right")
    ax7.set(yscale="log")
    ax7.set_ylim((10 ** 0, 10 ** 5))

    # cdr3 counts
    df["cdr3_count"] = df.groupby("cdr3")["cdr3"].transform("size")
    # x = df.groupby("cdr3").size().sort_values(ascending=False)
    ax8 = sns.countplot(
        x="cdr3_count",
        hue="gene",
        data=df.sort_values(by="cdr3_count", ascending=False),
        alpha=0.8,
        palette=pal,
    )
    ax8.set_ylabel("Count (log scale)")
    ax8.set_xlabel("Number of epitope partners per CDR3 sequence")
    ax8.set_title("CDR3 cross-reactivity")
    # ax8.set_xlabel('Number of occurrences of a given CDR3 in the dataset')
    # ax8.set_title('CDR3 promiscuity')
    ax8.legend(title="TCR chain", loc="upper right")
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
        for i in [ax5, ax6, ax7]
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


def assert_columns(df, input):  # noqa: A002
    assert df.columns.tolist() == [
        # "complex.id",
        # "Gene",
        # "CDR3",
        # "V",
        # "J",
        # "Species",
        # "MHC A",
        # "MHC B",
        # "MHC class",
        # "Epitope",
        # "Epitope gene",
        # "Epitope species",
        # "Reference",
        # "Method",
        # "Meta",
        # "CDR3fix",
        # "Score",
        "complex.id",
        "gene",
        "cdr3",
        "v.segm",
        "j.segm",
        "species",
        "mhc.a",
        "mhc.b",
        "mhc.class",
        "antigen.epitope",
        "antigen.gene",
        "antigen.species",
        "reference.id",
        "method",
        "meta",
        "cdr3fix",
        "vdjdb.score",
        "web.method",
        "web.method.seq",
        "web.cdr3fix.nc",
        "web.cdr3fix.unmp",
    ], f"Column headers of {input.resolve()} do not match the expected header names."


if __name__ == "__main__":

    # parse cli arguments
    args = create_parser()

    # create path objects for input and output files
    input_file = Path(args.input)
    output_file = Path(args.output)
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
    df = filter_vdjdb(
        input=input_file,
        tcr_chain=args.tcr_chain,
        species=args.species,
        drop_spurious=args.drop_spurious,
        mhc=args.mhc,
        hla=args.hla,
        specific_removal_epitope_reference=args.specific_removal_epitope_reference,
        specific_removal_references=args.specific_removal_references,
        keep_specific_references=args.keep_specific_references,
        keep_specific_epitopes=args.keep_specific_epitopes,
        length_restriction=args.length_restriction,
        downsample=args.downsample,
        terminal_only=args.terminal_only,
        middle_only=args.middle_only,
        terminal_replaced=args.terminal_replaced,
        middle_replaced=args.middle_replaced,
    )

    # create visualisation of dataset
    visualise_data(df, plot_file)

    # save entire dataset
    df.to_csv(
        (output_file.parent / (output_file.stem + "_full" + output_file.suffix)),
        index=False,
        sep=";",
    )
    logger.info(f"Saved processed full dataset to {output_file.resolve()}.")

    # save filtered output used for training
    # extract relevant columns
    columns = [
        "cdr3",
        "antigen.epitope",
    ]
    df = df.filter(items=columns)
    df.to_csv(output_file, index=False, sep=";")
    logger.info(f"Saved processed CDR3-epitope dataset to {output_file.resolve()}.")
