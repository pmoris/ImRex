import logging

import numpy as np
import pandas as pd


def sample_pairs(
    cdr3: str,
    df: pd.DataFrame,
    cdr3_column: str = "cdr3",
    epitope_column: str = "antigen.epitope",
) -> (str, str):
    """Sample an epitope for the given CDR3 sequence from the pool of other epitopes in the original positive dataset.

    NOTE: do not use a random_state for the sample function, since this will result in the same epitope
    being returned every time (for cdr3s with the same original epitope).

    Parameters
    ----------
    cdr3 : str
        The cdr3 sequence that should be matched with a negative epitope.
    df : pd.DataFrame
        A positive cdr3-epitope DataFrame with a "cdr3" and "antigen.epitope" column.
        Must have a class label column ("y") with "1" as the positive label.
    cdr3_column : str
        The header for the cdr3 column in the DataFrame.
    epitope_column : str
        The header for the epitope column in the DataFrame.

    Returns
    -------
    Tuple
        A tuple of a negative cdr3 and epitope sequence pair.
    """
    # check which epitopes occur as a positive partner for the current cdr3
    epitopes_to_exclude = df.loc[
        (df[cdr3_column] == cdr3) & (df["y"] == 1), epitope_column,
    ]
    # NOTE: for this to work, the original data source should either remain unmodified (to avoid epitopes paired with
    # the cdr3 as a negative example from showing up in this list), or by making sure the class labels are 1, in which
    # case the original dataframe should be given class labels before the sample_pairs function is called for the first time.

    # create pd.Series with all epitopes except for those that are positive partners of the current cdr3
    possible_epitopes = df.loc[
        ~df[epitope_column].isin(epitopes_to_exclude), epitope_column,
    ]

    # check if list is empty => cdr3 binds to every epitope present
    if possible_epitopes.empty:
        logger = logging.getLogger(__name__)
        logger.warning(
            f"CDR3 sequence {cdr3} is associated with every epitope in the dataset and will be discarded from the negatives."
        )
        return cdr3, np.NaN

    # sample 1 epitope from this list to pair with the cdr3 as a negative example
    else:
        sampled_epitope = possible_epitopes.sample(n=1).reset_index(drop=True)[0]
        return cdr3, sampled_epitope
