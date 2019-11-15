import os
import logging
import pandas as pd
import src.preprocessing.preprocess_vdjdb

# from io import StringIO
# import tempfile
from pathlib import Path

# file logger
log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logging.basicConfig(level=logging.INFO, format=log_fmt)
logger = logging.getLogger(__name__)


def test_filter_vdjdb():

    # csv_buffer = StringIO()
    # temp_test_file = tempfile.NamedTemporaryFile
    temp_test_file = "src/tests/preprocess_vdjdb_test.csv"

    # human TRB
    src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=Path("data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv"),
        output=Path(temp_test_file),
        tcr_chain="TRB",
        species="human",
        mhc="all",
        hla="all",
    )

    # df = pd.read_csv(csv_buffer, sep=";")
    df = pd.read_csv(temp_test_file, sep=";")
    assert df.shape == (32639, 2)
    os.remove(temp_test_file)

    # human TRA
    src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=Path("data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv"),
        output=Path(temp_test_file),
        tcr_chain="TRA",
        species="human",
        mhc="all",
        hla="all",
    )
    df = pd.read_csv(temp_test_file, sep=";")
    assert df.shape == (22248, 2)
    os.remove(temp_test_file)

    # all species all chains
    src.preprocessing.preprocess_vdjdb.filter_vdjdb(
        input=Path("data/raw/vdjdb/vdjdb-all-species-tra-trb-non-paired.tsv"),
        output=Path(temp_test_file),
        tcr_chain="all",
        species="all",
        mhc="all",
        hla="all",
    )
    df = pd.read_csv(temp_test_file, sep=";")
    assert df.shape == (59072, 2)
    os.remove(temp_test_file)
