import Levenshtein
import numpy as np
import pandas as pd


def calculate_distance(per_epitope_df: pd.DataFrame, train_df: pd.DataFrame):
    """
    calculate edit distances to training examples
    """
    epitope_dist_dict = {}

    for epitope in per_epitope_df["epitope"].unique():
        epitopes_to_check = train_df.loc[
            train_df["antigen.epitope"] != epitope, "antigen.epitope"
        ]
        distances_dict = {
            i: Levenshtein.distance(epitope, i) for i in epitopes_to_check.unique()
        }
        all_distances = epitopes_to_check.map(lambda x: distances_dict.get(x))

        epitope_dist_dict[epitope] = {
            "min_dist": min(distances_dict.values()),
            "median_dist": np.median(all_distances),
            "25th_quantile": np.percentile(all_distances, 25),
            "mean_dist": np.mean(all_distances),
        }

    per_epitope_df["min_dist"] = per_epitope_df["epitope"].map(
        lambda x: epitope_dist_dict.get(x).get("min_dist")
    )
    per_epitope_df["median_dist"] = per_epitope_df["epitope"].map(
        lambda x: epitope_dist_dict.get(x).get("median_dist")
    )
    per_epitope_df["25th_quantile"] = per_epitope_df["epitope"].map(
        lambda x: epitope_dist_dict.get(x).get("25th_quantile")
    )
    per_epitope_df["mean_dist"] = per_epitope_df["epitope"].map(
        lambda x: epitope_dist_dict.get(x).get("mean_dist")
    )

    return per_epitope_df
