""" Scenario for neural network with interaction map input. """
import argparse
import logging
from pathlib import Path

import Levenshtein
import numpy as np
import pandas as pd
import tensorflow as tf

from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.data.vdjdb_source import VdjdbSource
from src.scripts import io_helper
from src.neural import evaluation
from src.processing.data_stream import DataStream
from src.processing.padded_dataset_generator import padded_dataset_generator
from src.processing.separated_input_dataset_generator import (
    separated_input_dataset_generator,
)
from src.visualisation.plot import (
    derive_metrics_all,
    plot_all,
    roc_dist_corr,
    roc_per_epitope,
    roc_train_corr,
)


def create_parser():
    parser = argparse.ArgumentParser(
        description="Script to extract CDR3-epitope sequence pairs from VDJdb files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        dest="input",
        type=str,
        help="Directory containing directories named 'iteration_#' with h5 models and test datasets.",
    )
    # directory
    #   |- iteration 0
    #       |- metrics.csv
    #       |- roc.csv
    #       |- precision_recall.csv
    #       |- average_precision.csv
    #   |- iteration 1
    #       |- ...
    parser.add_argument(
        "--model_type",
        dest="model_type",
        choices=["padded", "separated"],
        type=str,
        help="Specifies the type of the model (padded interaction map or dual input type). Can be either 'padded' or 'separated'",
        required=True,
    )
    parser.add_argument(
        "--features",
        dest="features",
        type=str,
        help="A string of comma separated values listed in peptide_feature.featuresMap. Ignored for separate input models.",
        default="hydrophob,isoelectric,mass,hydrophil,charge",
    )
    parser.add_argument(
        "--operator",
        dest="operator",
        type=str,
        choices=["prod", "diff", "absdiff", "layer", "best"],
        help="Can be any of: prod, diff, absdiff, layer or best. Ignored for separate input models.",
        default="absdiff",
    )
    parser.add_argument(
        "--min_length_cdr3",
        dest="min_length_cdr3",
        type=int,
        help="Minimum CDR3 sequence length, used during negative reference filtering and padding.",
        default=None,
    )
    parser.add_argument(
        "--max_length_cdr3",
        dest="max_length_cdr3",
        type=int,
        help="Maximum CDR3 sequence length, used during negative reference filtering and padding.",
        default=None,
    )
    parser.add_argument(
        "--min_length_epitope",
        dest="min_length_epitope",
        type=int,
        help="Minimum epitope sequence length, used during negative reference filtering and padding.",
        default=None,
    )
    parser.add_argument(
        "--max_length_epitope",
        dest="max_length_epitope",
        type=int,
        help="Maximum epitope sequence length, used during negative reference filtering and padding.",
        default=None,
    )
    parser.add_argument(
        "--batch_size", dest="batch_size", type=int, help="Batch size.", default=128,
    )
    parser.add_argument(
        "--epitope_grouped",
        dest="epitope_grouped",
        action="store_true",
        help="Different plots should be generated for epitope-grouped models.",
        default=False,
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":

    # parse cli arguments
    args = create_parser()

    # create run name with timestamp
    run_name = io_helper.create_run_name("evaluate_test_folds")

    # check input directory
    input_directory = Path(args.input)
    assert input_directory.is_dir(), "Input is not a directory..."
    # create list of all "iterations_#" directories
    iterations_directories = [
        d
        for d in input_directory.iterdir()
        if d.is_dir() and d.name.startswith("iteration_")
    ]
    assert (
        iterations_directories
    ), "Could not find directories named 'iterations_#' in the input directory..."

    # create output directory
    output_dir = input_directory / run_name
    output_dir.mkdir(parents=True, exist_ok=False)

    # create logger in output directory
    io_helper.create_logger(run_name, evaluate_dir=output_dir)
    logger = logging.getLogger(__name__)

    # log arguments that were used
    for arg, value in sorted(vars(args).items()):
        logging.info("CLI argument %s: %r", arg, value)

    # create filepath and dataframe to store per epitope performance
    per_epitope_filepath = output_dir / "metrics_per_epitope.csv"
    per_epitope_df = pd.DataFrame()

    for iteration_dir in iterations_directories:

        # create metrics and figures (not required)
        # derive_metrics_all(iteration_dir, force=True)
        # plot_all(iteration_dir)

        # find latest model
        model_list = [model for model in iteration_dir.glob("*epoch*.h5")]
        model_path = sorted(
            model_list, key=lambda x: x.name.split("epoch")[1], reverse=True,
        )[0]

        # find test dataset
        data_path = list(iteration_dir.glob("test_fold_*.csv"))
        assert (
            len(data_path) == 1
        ), f"Found multiple files named 'test_fold_#.csv' inside {iteration_dir.name}..."

        # read test dataset
        data_source = VdjdbSource(
            filepath=data_path[0],
            headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
        )

        # filter on size
        min_length_cdr3 = (
            args.min_length_cdr3
            if args.min_length_cdr3
            else (data_source.data[data_source.headers["cdr3_header"]].str.len().min())
        )
        max_length_cdr3 = (
            args.max_length_cdr3
            if args.max_length_cdr3
            else (data_source.data[data_source.headers["cdr3_header"]].str.len().max())
        )
        min_length_epitope = (
            args.min_length_epitope
            if args.min_length_epitope
            else (
                data_source.data[data_source.headers["epitope_header"]].str.len().min()
            )
        )
        max_length_epitope = (
            args.max_length_epitope
            if args.max_length_epitope
            else (
                data_source.data[data_source.headers["epitope_header"]].str.len().max()
            )
        )

        # get list of features and operator based on input arguments
        if args.model_type == "padded":
            features_list = parse_features(args.features)
            operator = parse_operator(args.operator)
            feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)
            logger.info("features: " + str(features_list))
            logger.info("operator: " + str(operator))
        elif args.model_type == "separated":
            feature_builder = None

        # store range restrictions for cdr3 and epitope
        cdr3_range = (min_length_cdr3, max_length_cdr3)
        epitope_range = (min_length_epitope, max_length_epitope)
        logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
        logger.info(f"Filtered epitope sequences to length: {epitope_range}")

        if "y" not in data_source.data.columns:
            raise RuntimeError(
                "Dataset should contain positive and negative class labels."
            )

        # turn dataset into DataStream
        test_stream = DataStream(data_source)

        if args.model_type == "padded":
            dataset = padded_dataset_generator(
                data_stream=test_stream,
                feature_builder=feature_builder,
                cdr3_range=cdr3_range,
                epitope_range=epitope_range,
                neg_shuffle=False,
                export_path=None,
            )
        elif args.model_type == "separated":
            dataset = separated_input_dataset_generator(
                data_stream=test_stream,
                cdr3_range=cdr3_range,
                epitope_range=epitope_range,
                neg_shuffle=False,
                export_path=None,
            )
        dataset = dataset.batch(args.batch_size)

        # load model
        model = tf.keras.models.load_model(model_path)

        # evaluate
        try:
            metrics_df = evaluation.evaluate_model(model=model, dataset=dataset)
        except ValueError as e:
            raise ValueError(
                "Make sure the correct model type and padding lengths are specified. The latter can be omitted if every fold contains an example of the shortest and longest length, otherwise they should be provided as input arguments."
            ) from e
        metrics_filepath = output_dir / f"metrics_{iteration_dir.name}.csv"
        metrics_df.to_csv(metrics_filepath, index=False)
        logger.info(
            f"Saved metrics for {iteration_dir.name} in {metrics_filepath.absolute()}."
        )

        predictions_df = evaluation.predictions_model(model=model, dataset=dataset)
        predictions_filepath = output_dir / f"predictions_{iteration_dir.name}.csv"
        predictions_df.to_csv(predictions_filepath, index=False)
        logger.info(
            f"Saved predictions for {iteration_dir.name} in {predictions_filepath.absolute()}."
        )

        try:
            iteration_metrics_df = evaluation.evaluate_per_epitope(
                model=model,
                data_source=data_source,
                batch_size=args.batch_size,
                model_type=args.model_type,
                feature_builder=feature_builder,
                cdr3_range=cdr3_range,
                epitope_range=epitope_range,
            )
        except ValueError as e:
            raise ValueError(
                "Make sure the correct model type and padding lengths are specified. The latter can be omitted if every fold contains an example of the shortest and longest length, otherwise they should be provided as input arguments."
            ) from e
        iteration_metrics_df["iteration"] = iteration_dir.name

        # add training dataset size per epitope
        train_path = list(iteration_dir.glob("train_fold_*.csv"))
        assert (
            len(data_path) == 1
        ), f"Found multiple files named 'train_fold_#.csv' inside {iteration_dir.name}..."
        train_df = pd.read_csv(train_path[0], sep=";")
        iteration_metrics_df["train_size"] = iteration_metrics_df.apply(
            lambda x: np.sum(train_df["antigen.epitope"] == x["epitope"]), axis=1
        )

        # calculate edit distances to training examples
        epitope_dist_dict = {}
        for epitope in iteration_metrics_df["epitope"].unique():
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

        iteration_metrics_df["min_dist"] = iteration_metrics_df["epitope"].map(
            lambda x: epitope_dist_dict.get(x).get("min_dist")
        )
        iteration_metrics_df["median_dist"] = iteration_metrics_df["epitope"].map(
            lambda x: epitope_dist_dict.get(x).get("median_dist")
        )
        iteration_metrics_df["25th_quantile"] = iteration_metrics_df["epitope"].map(
            lambda x: epitope_dist_dict.get(x).get("25th_quantile")
        )
        iteration_metrics_df["mean_dist"] = iteration_metrics_df["epitope"].map(
            lambda x: epitope_dist_dict.get(x).get("mean_dist")
        )

        # combine df for current iteration with list of all iterations
        per_epitope_df = per_epitope_df.append(iteration_metrics_df)

    per_epitope_df.to_csv(per_epitope_filepath, index=False)
    logger.info(f"Saved per-epitope metrics in {per_epitope_filepath.absolute()}.")

    # create plots
    roc_per_epitope(
        eval_df=per_epitope_df,
        output_path=output_dir / "roc_per_epitope.pdf",
        min_obs=30,
        min_iterations=5,
        grouped=args.epitope_grouped,
    )

    roc_train_corr(
        eval_df=per_epitope_df, output_path=output_dir / "roc_train_correlation.pdf"
    )

    roc_dist_corr(
        eval_df=per_epitope_df, output_path=output_dir / "roc_dist_correlation.pdf"
    )
