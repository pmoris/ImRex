import logging

import pandas as pd


def export_data(x_seq, y_seq, export_path):
    """Export dataset to csv. """
    df = pd.DataFrame(x_seq, columns=["cdr3", "antigen.epitope"])
    df["class"] = y_seq

    logger = logging.getLogger(__name__)
    logger.info(f"Saving train/test fold in: {export_path}")
    df.to_csv(export_path, sep=";", index=False)
