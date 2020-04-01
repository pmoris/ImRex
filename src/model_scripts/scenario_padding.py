""" Scenario for neural network with interaction map input. """
import argparse
import logging

from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.control_cdr3_source import ControlCDR3Source
from src.data.vdjdb_source import VdjdbSource
from src.model_scripts import pipeline
from src.neural.trainer import get_output_path, Trainer
from src.processing.cv_folds import cv_splitter
from src.processing.inverse_map import InverseMap
from src.processing.padded_dataset_generator import padded_dataset_generator
from src.processing.splitter import splitter


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
        "--batch_size", dest="batch_size", type=int, help="Batch size.", default=128,
    )
    parser.add_argument(
        "--epochs",
        dest="epochs",
        type=int,
        help="Number of epochs to train for.",
        default=40,
    )
    parser.add_argument(
        "--neg_ref",
        dest="neg_ref",
        type=str,
        help="Whether to generate negatives from CDR3 reference sequences or by shuffling positive examples.",
        default=None,
    )
    parser.add_argument(
        "--neg_ratio",
        dest="neg_ratio",
        type=float,
        help="Proportion of positive to negative samples.",
        default=0.5,
    )
    parser.add_argument(
        "--val_split",
        dest="val_split",
        type=float,
        help="The proportion of the dataset to include in the test split. Omit to use cross-validation.",
        default=None,
    )
    parser.add_argument(
        "--epitope_grouped_cv",
        dest="epitope_grouped_cv",
        action="store_true",
        help="When val_split is omitted, indicates whether to use an epitope-grouped cv over a normal k-fold cv.",
        default=False,
    )
    parser.add_argument(
        "--one_epitope_out_cv",
        dest="one_epitope_out_cv",
        action="store_true",
        help="When val_split is omitted and epitope_grouped_cv is selected, indicates whether to use leave-1-epitope-out cv",
        default=False,
    )
    parser.add_argument(
        "--neg_gen_full",
        dest="neg_gen_full",
        action="store_true",
        help="When selected, the entire dataset will be shuffled to generate negatives, which can cause individual cdr3/epitope sequences to reappear in both train and test sets. Otherwise, negatives will be generated through shuffling within each train/test set. Otherwise,",
        default=False,
    )
    parser.add_argument(
        "--n_folds",
        dest="n_folds",
        type=int,
        help="Number of folds to use during cross-validation.",
        default=5,
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
        "--name",
        dest="name",
        type=str,
        help="Name under which the model and log files will be stored, appended with the date and time.",
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
        "--early_stop",
        dest="early_stop",
        action="store_true",
        help="Terminate model training when the performance metric stops improving (uses tf.keras.callbacks.EarlyStopping).",
        default=False,
    )
    parser.add_argument(
        "--include_learning_rate_reduction",
        dest="include_learning_rate_reduction",
        action="store_true",
        help="Reduce the learning rate when the performance metric has stopped improving (tf.keras.callbacks.ReduceLROnPlateau).",
        default=False,
    )
    parser.add_argument(
        "--optimizer",
        dest="optimizer",
        type=str,
        choices=["rmsprop", "adam", "SGD"],
        help="Can be any of: rmsprop, adam or SGD.",
        default="rmsprop",
    )
    parser.add_argument(
        "--learning_rate",
        dest="learning_rate",
        type=float,
        help="Learning rate supplied to the selected optimizer",
        default=None,
    )
    parser.add_argument(
        "--model",
        dest="model",
        type=str,
        choices=["model_padded", "model_selu", "model_small", "model_small_selu"],
        help="The type of padded model to use.",
        default="model_padded",
    )
    parser.add_argument(
        "--depth1",
        dest="depth1",
        type=int,
        help="Depth of the first convolutional layer. Only used for the small models.",
    )
    parser.add_argument(
        "--depth2",
        dest="depth2",
        type=int,
        help="Depth of the second convolutional layer. Only used for the small models.",
    )
    parser.add_argument(
        "--dropout1",
        dest="dropout1",
        type=float,
        help="Dropout after the first convolutional layer. Only used for the small models.",
    )
    parser.add_argument(
        "--dropout2",
        dest="dropout2",
        type=float,
        help="Dropout after the second convolutional layer. Only used for the small models.",
    )
    parser.add_argument(
        "--gap",
        dest="gap",
        action="store_true",
        help="Use global average max pool. Only used for the small models.",
        default=False,
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    # parse cli arguments
    args = create_parser()

    # create logger and log file
    run_name = pipeline.create_run_name(args.name)
    pipeline.create_logger(run_name)
    logger = logging.getLogger(__name__)

    # log arguments that were used
    for arg, value in sorted(vars(args).items()):
        logging.info("CLI argument %s: %r", arg, value)

    # read (positive) data
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

    # check argument compatability
    if args.epitope_grouped_cv and args.val_split is not None:
        raise RuntimeError("Cannot test epitope-grouped without k folds.")

    logger.info("features: " + str(features_list))
    logger.info("operator: " + str(operator))
    logger.info("neg_ref: " + str(args.neg_ref))
    logger.info("epitope_grouped_cv: " + str(args.epitope_grouped_cv))
    logger.info("neg_gen_full: " + str(args.neg_gen_full))

    inverse_map = InverseMap()

    # store range restrictions for cdr3 and epitope
    cdr3_range = (args.min_length_cdr3, args.max_length_cdr3)
    epitope_range = (args.min_length_epitope, args.max_length_epitope)
    logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
    logger.info(f"Filtered epitope sequences to length: {epitope_range}")

    trainer = Trainer(
        args.epochs,
        include_learning_rate_reduction=args.include_learning_rate_reduction,
        include_early_stop=args.early_stop,
        lookup=inverse_map,
    )

    if args.model == "model_padded":
        from src.models.model_padded import ModelPadded

        model = ModelPadded(
            width=args.max_length_cdr3,
            height=args.max_length_epitope,
            name=run_name,
            channels=feature_builder.get_number_layers(),
            optimizer=args.optimizer,
            learning_rate=args.learning_rate,
        )
    elif args.model == "model_selu":
        from src.models.model_padded_selu import ModelPaddedSelu

        model = ModelPaddedSelu(
            width=args.max_length_cdr3,
            height=args.max_length_epitope,
            name=run_name,
            channels=feature_builder.get_number_layers(),
            optimizer=args.optimizer,
            learning_rate=args.learning_rate,
        )
    elif args.model == "model_small":
        from src.models.model_padded_small import ModelPaddedSmall

        model = ModelPaddedSmall(
            width=args.max_length_cdr3,
            height=args.max_length_epitope,
            name=run_name,
            channels=feature_builder.get_number_layers(),
            optimizer=args.optimizer,
            learning_rate=args.learning_rate,
            depth1=args.depth1,
            depth2=args.depth2,
            dropout1=args.dropout1,
            dropout2=args.dropout2,
            gap=args.gap,
        )
    elif args.model == "model_small_selu":
        from src.models.model_padded_small_selu import ModelPaddedSmallSelu

        model = ModelPaddedSmallSelu(
            width=args.max_length_cdr3,
            height=args.max_length_epitope,
            name=run_name,
            channels=feature_builder.get_number_layers(),
            optimizer=args.optimizer,
            learning_rate=args.learning_rate,
            depth1=args.depth1,
            depth2=args.depth2,
            dropout1=args.dropout1,
            dropout2=args.dropout2,
            gap=args.gap,
        )

    logger.info(f"Built model {model.base_name}:")
    # model.summary() is logged inside trainer.py

    # if negative reference dataset is provided, draw negatives from it
    if args.neg_ref:
        logger.info(f"Generating negative examples from negative reference CDR3 set.")
        negative_source = ControlCDR3Source(
            filepath=args.neg_ref,
            min_length=args.min_length_cdr3,
            max_length=args.max_length_cdr3,
        )
        data_source.generate_negatives_from_ref(negative_source)
        neg_shuffle = False

    # otherwise, generate negatives through shuffling
    else:
        # if neg_gen_full is True, generate negatives once on the entire dataset
        if args.neg_gen_full:
            logger.info(
                f"Generating negative examples through shuffling on the entire dataset prior to train/test fold creation."
            )
            data_source.generate_negatives()
            neg_shuffle = False

        # otherwise, generate negatives within each train/test set during tf dataset creation
        else:
            logger.info(
                f"Generating negative examples through shuffling within each train/test fold."
            )
            neg_shuffle = True

    # if a fixed train-test split ratio is provided...
    if args.val_split is not None:
        train, val = splitter(data_source, test_size=args.val_split)
        iterations = [(train, val)]

    # ...otherwise use the specified cross validation scheme
    else:
        iterations = cv_splitter(
            data_source=data_source,
            n_folds=args.n_folds,
            epitope_grouped=args.epitope_grouped_cv,
            one_out=args.one_epitope_out_cv,
        )

    for iteration, (train, val) in enumerate(iterations):

        logger.info(f"Iteration: {iteration}")
        logger.info(f"Batch size: {args.batch_size}")
        logger.info(f"Train set: {len(train)}")
        logger.info(f"Test set: {len(val)}")

        # retrieve model output directory and create data directory to store generated datasets with positive and negative examples
        train_fold_output, test_fold_output = (
            get_output_path(
                base_name=run_name,
                file_name=f"train_fold_{iteration}.csv",
                iteration=iteration,
            ),
            get_output_path(
                base_name=run_name,
                file_name=f"test_fold_{iteration}.csv",
                iteration=iteration,
            ),
        )

        train_data = padded_dataset_generator(
            data_stream=train,
            feature_builder=feature_builder,
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
            neg_shuffle=neg_shuffle,
            export_path=train_fold_output,
        )
        val_data = padded_dataset_generator(
            data_stream=val,
            feature_builder=feature_builder,
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
            inverse_map=inverse_map,
            neg_shuffle=neg_shuffle,
            export_path=test_fold_output,
        )

        # get length of train dataset
        train_length = len(train) if not neg_shuffle else len(train) * 2

        # shuffle and batch train data
        train_data = train_data.shuffle(
            # buffer equals size of dataset, because positives and negatives are grouped
            buffer_size=train_length,
            seed=42,
            # reshuffle to make each epoch see a different order of examples
            reshuffle_each_iteration=True,
        ).batch(args.batch_size)

        # batch validation data
        val_data = val_data.batch(args.batch_size)

        trainer.train(model, train_data, val_data, iteration=iteration)
