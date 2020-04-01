""" Scenario for neural network with interaction map input. """
import argparse
import logging

import tensorflow as tf

from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.control_cdr3_source import ControlCDR3Source
from src.data.vdjdb_source import VdjdbSource
from src.model_scripts import pipeline
from src.neural import evaluation
from src.processing.data_stream import DataStream
from src.processing.padded_dataset_generator import padded_dataset_generator


def create_parser():
    parser = argparse.ArgumentParser(
        description="Script to extract CDR3-epitope sequence pairs from VDJdb files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--data_path",
        dest="data_path",
        type=str,
        help="Input csv dataset, as supplied by preprocess_vdjdb script.",
        default=PROJECT_ROOT
        / "data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x.csv",
    )
    parser.add_argument(
        "--model",
        dest="model",
        type=str,
        help="Path to the trained model in h5 format.",
        required=True,
    )
    parser.add_argument(
        "--model_type",
        dest="model_type",
        choices=["padded", "separated"],
        type=str,
        help="Specifies the type of the model (padded interaction map or dual input type). Can be either 'padded' or 'separated'",
        required=True,
    )
    parser.add_argument(
        "--neg_ref",
        dest="neg_ref",
        type=str,
        help="Generate negatives from CDR3 reference sequences or by shuffling positive examples.",
        default=None,
    )
    parser.add_argument(
        "--neg_gen",
        dest="neg_gen",
        action="store_true",
        help="Generate negatives through shuffling the positive pairs in the supplied dataset.",
        default=False,
    )
    # parser.add_argument(
    #     "--neg_augment",
    #     dest="neg_augment",
    #     type=str,
    #     help="The path to a negative reference cdr3 set, used to augment the shuffled negatives with additional examples. Ignored by default. Must be used in conjunction with the --augment_amount argument",
    #     default=None,
    # )
    # parser.add_argument(
    #     "--augment_amount",
    #     dest="augment_amount",
    #     type=int,
    #     help="The number of additional negatives to generate from the negative reference cdr3 set. Must be used in conjunction with the --neg_augment argument.",
    #     default=None,
    # )
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
        "--name",
        dest="name",
        type=str,
        help="Name under which the test output and log files will be stored, appended with the date and time.",
        required=True,
    )
    parser.add_argument(
        "--features",
        dest="features",
        type=str,
        help="A string of comma separated values listed in peptide_feature.featuresMap.",
        default="hydrophob,isoelectric,mass,hydrophil,charge",
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
        "--batch_size", dest="batch_size", type=int, help="Batch size.", default=128,
    )
    parser.add_argument(
        "--per_epitope",
        dest="per_epitope",
        action="store_true",
        help="Perform all evaluations on a per-epitope basis.",
        default=False,
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":

    # parse cli arguments
    args = create_parser()

    # create logger and log file
    run_name = pipeline.create_run_name(args.name)
    output_dir = pipeline.create_evaluate_path(args.name, args.model)
    pipeline.create_logger(run_name, evaluate_dir=output_dir)
    logger = logging.getLogger(__name__)

    # log arguments that were used
    for arg, value in sorted(vars(args).items()):
        logging.info("CLI argument %s: %r", arg, value)

    # read positive and negative data
    data_source = VdjdbSource(
        filepath=args.data_path,
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )

    # filter on size
    data_source.length_filter(
        args.min_length_cdr3,
        args.max_length_cdr3,
        args.min_length_epitope,
        args.max_length_epitope,
    )

    # get list of features and operator based on input arguments
    features_list = parse_features(args.features)
    operator = parse_operator(args.operator)
    feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)

    logger.info("features: " + str(features_list))
    logger.info("operator: " + str(operator))
    logger.info("neg_ref: " + str(args.neg_ref))
    logger.info("neg_gen: " + str(args.neg_gen))

    # store range restrictions for cdr3 and epitope
    cdr3_range = (args.min_length_cdr3, args.max_length_cdr3)
    epitope_range = (args.min_length_epitope, args.max_length_epitope)
    logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
    logger.info(f"Filtered epitope sequences to length: {epitope_range}")

    # if negative reference dataset is provided, draw negatives from it
    if args.neg_ref:
        logger.info(f"Generating negative examples from negative reference CDR3 set.")
        if "y" in data_source.data.columns:
            if data_source.data.loc[data_source.data["y"] == 0].shape[0] > 0:
                raise RuntimeError(
                    "Dataset already contains negative class labels. Do not use --neg_ref argument."
                )
        # read negative reference data
        negative_source = ControlCDR3Source(
            filepath=args.neg_ref,
            min_length=args.min_length_cdr3,
            max_length=args.max_length_cdr3,
        )
        data_source.add_pos_labels()
        data_source.generate_negatives_from_ref(negative_source)

    # or generate negatives through shuffling
    elif args.neg_gen:
        logger.info(
            f"Generating negative examples through shuffling the positive dataset."
        )
        if "y" in data_source.data.columns:
            if data_source.data.loc[data_source.data["y"] == 0].shape[0] > 0:
                raise RuntimeError(
                    "Dataset already contains negative class labels. Do not use --neg_gen argument."
                )
        data_source.add_pos_labels()
        data_source.generate_negatives()

    # otherwise negatives should be present in dataset already
    else:
        if "y" not in data_source.data.columns:
            raise RuntimeError(
                "Dataset should contain positive and negative class labels if neither --neg_ref or --neg_gen arguments are used."
            )
        logger.info(
            f"Not generating negatives. The supplied dataset already contains both positives and negatives: {data_source.data.groupby('y').size()}"
        )

    # turn dataset into DataStream
    test_stream = DataStream(data_source)

    # retrieve evaluation output directory and create filepath store generated datasets
    evaluation_output_dataset = output_dir / "evaluation_dataset.csv"

    val_data = padded_dataset_generator(
        data_stream=test_stream,
        feature_builder=feature_builder,
        cdr3_range=cdr3_range,
        epitope_range=epitope_range,
        neg_shuffle=False,
        export_path=evaluation_output_dataset,
    )
    val_data = val_data.batch(args.batch_size)

    # load model
    model = tf.keras.models.load_model(args.model)

    # evaluate
    evaluation.evaluate_model(model=model, dataset=val_data, output_dir=output_dir)
    evaluation.predictions_model(model=model, dataset=val_data, output_dir=output_dir)

    if args.per_epitope:
        evaluation.evaluate_per_epitope(
            model=model,
            data_source=data_source,
            output_dir=output_dir,
            batch_size=args.batch_size,
            model_type=args.model_type,
            feature_builder=feature_builder,
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
        )
