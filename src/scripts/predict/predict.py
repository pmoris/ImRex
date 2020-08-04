import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import tensorflow as tf

from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.neural import evaluation
from src.processing.data_stream import DataStream
from src.processing.image_generator import ImageGenerator
from src.processing.image_padding import ImagePadding
from src.processing.padded_dataset_generator import padded_dataset_generator
from src.processing.separated_input_dataset_generator import (
    BlosumImageGenerator,
    BlosumPadding,
    separated_input_dataset_generator,
)


def create_parser():
    parser = argparse.ArgumentParser(
        description="Script to predict binding between CDR3 and epitope sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        dest="input",
        type=str,
        help="Input dataset with paired CDR3 and epitope sequences, in .csv format. Default column labels are 'cdr3' and 'antigen.epitope' and a ';' separator.",
        required=True,
    )
    parser.add_argument(
        "--separator",
        dest="sep",
        type=str,
        help="Separator (delimiter) used for the input dataset. Default is ';'.",
        default=";",
    )
    parser.add_argument(
        "--cdr3_column",
        dest="cdr3_column",
        type=str,
        help="CDR3 column label, default is 'cdr3'.",
        default="cdr3",
    )
    parser.add_argument(
        "--epitope_column",
        dest="epitope_column",
        type=str,
        help="CDR3 column label, default is 'antigen.epitope'.",
        default="antigen.epitope",
    )
    parser.add_argument(
        "--output",
        dest="output",
        type=str,
        help="Specify the output file.",
        required=True,
    )
    parser.add_argument(
        "--model",
        dest="model",
        type=str,
        help="Path to the trained model in h5 format, see ./models/models-pretrained. Depending on the selected model the other parameters (such as min/max sequence length) will need to be adjusted.",
        default=PROJECT_ROOT
        / "models/models-pretrained/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001/2020-07-24_19-18-39_trbmhcidown-shuffle-padded-b32-lre4-reg001.h5",
    )
    parser.add_argument(
        "--model_type",
        dest="model_type",
        choices=["padded", "separated"],
        type=str,
        help="Specifies the type of the model (padded interaction map or dual input type). Can be either 'padded' (= interaction map, default) or 'separated'",
        default="padded",
    )
    parser.add_argument(
        "--min_length_cdr3",
        dest="min_length_cdr3",
        type=int,
        help="Minimum CDR3 sequence length, used during negative reference filtering and padding.",
        default=10,
    )
    parser.add_argument(
        "--max_length_cdr3",
        dest="max_length_cdr3",
        type=int,
        help="Maximum CDR3 sequence length, used during negative reference filtering and padding.",
        default=20,
    )
    parser.add_argument(
        "--min_length_epitope",
        dest="min_length_epitope",
        type=int,
        help="Minimum epitope sequence length, used during negative reference filtering and padding.",
        default=8,
    )
    parser.add_argument(
        "--max_length_epitope",
        dest="max_length_epitope",
        type=int,
        help="Maximum epitope sequence length, used during negative reference filtering and padding.",
        default=11,
    )
    parser.add_argument(
        "--features",
        dest="features",
        type=str,
        help="A string of comma separated values listed in peptide_feature.featuresMap. Used for padded (interaction map) type models.",
        default="hydrophob,isoelectric,mass,hydrophil",
    )
    parser.add_argument(
        "--operator",
        dest="operator",
        type=str,
        choices=["prod", "diff", "absdiff", "layer", "best"],
        help="Can be any of: prod, diff, absdiff, layer or best.",
        default="absdiff",
    )
    parser.add_argument(
        "--batch_size",
        dest="batch_size",
        type=int,
        help="Batch size used during prection, can be increased depending on GPU memory.",
        default=128,
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # parse cli arguments
    args = create_parser()

    output = Path(args.output)
    assert output.parent.exists(), "Provided output directory does not exist."

    log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    # suppress tf logging
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"  # ERROR
    logging.getLogger("tensorflow").setLevel(logging.ERROR)
    logger = logging.getLogger(__name__)

    # log arguments that were used
    for arg, value in sorted(vars(args).items()):
        logging.info("CLI argument %s: %r", arg, value)

    # read in dataset
    data_source = VdjdbSource(
        filepath=args.input,
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
        sep=args.sep,
    )

    # filter on size
    data_source.length_filter(
        args.min_length_cdr3,
        args.max_length_cdr3,
        args.min_length_epitope,
        args.max_length_epitope,
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
    cdr3_range = (args.min_length_cdr3, args.max_length_cdr3)
    epitope_range = (args.min_length_epitope, args.max_length_epitope)
    logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
    logger.info(f"Filtered epitope sequences to length: {epitope_range}")

    data_stream = DataStream(
        data_source.data[
            [data_source.headers["cdr3_header"], data_source.headers["epitope_header"]]
        ].to_numpy()
    )

    # store maximum cdr3 and epitope sequence length as width and height of the 2d arrays
    cdr3_max = cdr3_range[1]
    epitope_max = epitope_range[1]

    # setup tf datasets
    if args.model_type == "padded":

        image_gen = ImageGenerator(data_stream, feature_builder, has_label=False)
        image_padding = ImagePadding(
            image_gen, cdr3_max, epitope_max, pad_value=0, has_label=False
        )

        padding_array = np.array([image for image in image_padding])

        dataset = tf.data.Dataset.from_tensor_slices(padding_array)
        dataset = dataset.batch(args.batch_size)

    elif args.model_type == "separated":
        # create BLOSUM encoded arrays
        blosum_encoding = BlosumImageGenerator(data_stream, has_label=False)

        # pad with zero-columns
        blosum_padding = BlosumPadding(
            blosum_encoding, cdr3_max, epitope_max, pad_value=0, has_label=False
        )

        cdr3, epitope = zip(*blosum_padding)

        dataset = (np.array(cdr3), np.array(epitope))

    # load model
    model = tf.keras.models.load_model(args.model)

    # predictions
    try:
        y_pred = model.predict(dataset)
        data_source.data["prediction_score"] = y_pred
    except ValueError as e:
        raise ValueError(
            "Make sure the correct model type, optional features and padding lengths are specified."
        ) from e

    # save output
    data_source.data.to_csv(output, index=False)
    logger.info(f"Saved predictions in {output.absolute()}.")
