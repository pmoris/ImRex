# Test to see if scripts and datasets result in the expected number of entries
# Import notes:
# Duplicates are always removed as the very last step, based only on the cdr3 and epitope columns.
# When not specified explicitly, spurious entries are removed by default.

# import os
# import logging

# from pathlib import Path
# import pandas as pd
# from io import StringIO
# import tempfile

from src.config import PROJECT_ROOT
import src.preprocessing.preprocess_vdjdb

# enable file logger because preprocessing function expects it...
# log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
# logging.basicConfig(level=logging.INFO, format=log_fmt)
# logger = logging.getLogger(__name__)
# disable logger from showing up on console
# logger.disabled = True


def test_filter_vdjdb():

    #####################
    # vdjdb web release #
    #####################

    # NO LONGER FUNCTIONAL DUE TO CHANGES IN COLUMN NAMES
    # deprecated code is kept as log on how to use temporary testing files

    # # csv_buffer = StringIO()
    # # temp_test_file = tempfile.NamedTemporaryFile
    # # temp_test_file = "src/tests/preprocess_vdjdb_test.csv"

    # # human TRB - web
    # df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
    #     # input=Path("data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv"),
    #     input=Path("data/raw/vdjdb/vdjdb-browser.tsv"),
    #     # output=Path(temp_test_file),
    #     tcr_chain="TRB",
    #     species="human",
    #     mhc="all",
    #     hla="all",
    # )
    # # df = pd.read_csv(csv_buffer, sep=";")
    # # df = pd.read_csv(temp_test_file, sep=";")
    # assert df.shape == (32639, 2)
    # # os.remove(temp_test_file)

    # # human TRA - web
    # df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
    #     input=Path("data/raw/vdjdb/vdjdb-browser.tsv"),
    #     tcr_chain="TRA",
    #     species="human",
    #     mhc="all",
    #     hla="all",
    # )
    # assert df.shape == (22248, 2)

    # # all species all chains - web
    # df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
    #     input=Path("data/raw/vdjdb/vdjdb-browser.tsv"),
    #     tcr_chain="all",
    #     species="all",
    #     mhc="all",
    #     hla="all",
    # )
    # assert df.shape == (59072, 2)

    ########################
    # vdjdb github release #
    ########################

    # human TRB with spurious - vdjdb
    df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=PROJECT_ROOT / "data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt",
        tcr_chain="TRB",
        species="human",
        mhc="all",
        hla="all",
    )
    assert df.shape == (32639, 2)

    # human TRB without spurious - vdjdb
    df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=PROJECT_ROOT / "data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt",
        tcr_chain="TRB",
        species="human",
        mhc="all",
        hla="all",
        drop_spurious=False,
    )
    assert df.shape == (33028, 2)

    # human TRB without spurious without 10x - vdjdb
    df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=PROJECT_ROOT / "data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt",
        tcr_chain="TRB",
        species="human",
        mhc="all",
        hla="all",
        specific_removal_epitope_reference=["KLGGALQAK", "10xgenomics"],
    )
    assert df.shape == (20024, 2)

    # human TRA - vdjdb
    df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=PROJECT_ROOT / "data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt",
        tcr_chain="TRA",
        species="human",
        mhc="all",
        hla="all",
        drop_spurious=False,
    )
    assert df.shape == (23564, 2)

    # all species all chains with spurious - vdjdb
    df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=PROJECT_ROOT / "data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt",
        tcr_chain="all",
        species="all",
        mhc="all",
        hla="all",
        drop_spurious=False,
    )
    assert df.shape == (61047, 2)

    # all species all chains - vdjdb
    df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=PROJECT_ROOT / "data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt",
        tcr_chain="all",
        species="all",
        mhc="all",
        hla="all",
    )
    assert df.shape == (59072, 2)

    # human TRB length 10-20cdr3 10-13epitope - vdjdb - no10x - no spurious
    df = src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=PROJECT_ROOT / "data/raw/vdjdb/vdjdb-2019-08-08/vdjdb.txt",
        tcr_chain="all",
        species="human",
        drop_spurious=True,
        mhc="all",
        hla="all",
        specific_removal_epitope_reference=None,
        specific_removal_references=["10xgenomics"],
        keep_specific_references=None,
        length_restriction=["10", "20", "8", "13"],
    )
    assert df.shape == (20646, 2)
