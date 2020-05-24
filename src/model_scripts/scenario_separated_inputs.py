""" Scenario for neural network with dual sequence input. """
import argparse
import logging

from src.config import PROJECT_ROOT
from src.data.control_cdr3_source import ControlCDR3Source
from src.data.vdjdb_source import VdjdbSource
from src.model_scripts import pipeline
from src.neural.trainer import get_output_path, Trainer
from src.processing.cv_folds import cv_splitter
from src.processing.separated_input_dataset_generator import (
    separated_input_dataset_generator,
)
from src.processing.splitter import splitter


def create_parser():
    parser = argparse.ArgumentParser(
        description="Script to train a separated inputs cdr3-epitpe prediction model.",
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
        "--cv",
        dest="cross_validation",
        type=str,
        help="The type of cross-validation strategy to use. Can be any of: ",
        choices=[
            "kfold",
            "repeatedkfold",
            "epitope_grouped",
            "one_epitope_out",
            "epitope_grouped_shuffle",
        ],
        default=None,
    )
    parser.add_argument(
        "--neg_gen_full",
        dest="neg_gen_full",
        action="store_true",
        help="When selected, the entire dataset will be shuffled to generate negatives, which can cause individual cdr3/epitope sequences to reappear in both train and test sets. Otherwise, negatives will be generated through shuffling within each train/test set.",
        default=False,
    )
    parser.add_argument(
        "--neg_augment",
        dest="neg_augment",
        type=str,
        help="The path to a negative reference cdr3 set, used to augment the shuffled negatives with additional examples. Ignored by default. Must be used in conjunction with the --augment_amount argument",
        default=None,
    )
    parser.add_argument(
        "--augment_amount",
        dest="augment_amount",
        type=int,
        help="The number of additional negatives to generate from the negative reference cdr3 set. Must be used in conjunction with the --neg_augment argument.",
        default=None,
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
        choices=["separated", "nettcr", "nettcr_selu",],
        help="The type of separated input model to use.",
        default="separated",
    )
    parser.add_argument(
        "--regularization",
        dest="regularization",
        type=float,
        help="Regularization applied to all layers.",
        default=None,
    )
    parser.add_argument(
        "--dropout_conv",
        dest="dropout_conv",
        type=float,
        help="Dropout rate applied applied to the convolutional layers.",
        default=None,
    )
    parser.add_argument(
        "--dropout_dense",
        dest="dropout_dense",
        type=float,
        help="Dropout rate applied applied to the dense layers.",
        default=None,
    )
    parser.add_argument(
        "--activation_function",
        dest="activation_function",
        type=str,
        help="The activation function used for all layers except for the output.",
        default="selu",
    )
    parser.add_argument(
        "--disable_file_log",
        dest="disable_file_log",
        action="store_false",
        help="Disable logging to file.",
        default=True,
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    # parse cli arguments
    args = create_parser()

    # create logger and log file
    run_name = pipeline.create_run_name(args.name)
    pipeline.create_logger(run_name, log_to_file=args.disable_file_log)
    logger = logging.getLogger(__name__)

    # log arguments that were used
    for arg, value in sorted(vars(args).items()):
        logging.info("CLI argument %s: %r", arg, value)

    # read (positive) data
    data_source = VdjdbSource(
        filepath=args.data_path,
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )
    data_source.add_pos_labels()

    # filter on size
    data_source.length_filter(
        args.min_length_cdr3,
        args.max_length_cdr3,
        args.min_length_epitope,
        args.max_length_epitope,
    )

    # check argument compatability
    if args.val_split and args.cross_validation:
        raise RuntimeError(
            "Invalid arguments. Val_split should not be used in conjunction with a cv strategy."
        )
    elif bool(args.neg_augment) ^ bool(args.augment_amount):
        raise RuntimeError(
            "If negatives should be augmented, both the neg_agument and augment_amount should be supplied."
        )
    elif bool(args.n_folds) ^ bool(args.cross_validation):
        raise RuntimeError(
            "For cross-validation, both the number of folds and the type of cv should be specified."
        )

    logger.info("neg_ref: " + str(args.neg_ref))
    logger.info("cv_type: " + str(args.cross_validation))
    logger.info("neg_gen_full: " + str(args.neg_gen_full))

    # store range restrictions for cdr3 and epitope
    cdr3_range = (args.min_length_cdr3, args.max_length_cdr3)
    epitope_range = (args.min_length_epitope, args.max_length_epitope)
    logger.info(f"Filtered CDR3 sequences to length: {cdr3_range}")
    logger.info(f"Filtered epitope sequences to length: {epitope_range}")

    trainer = Trainer(
        args.epochs,
        include_learning_rate_reduction=args.include_learning_rate_reduction,
        include_early_stop=args.early_stop,
        verbose=args.disable_file_log,
    )

    if args.model == "separated":
        from src.models.model_separated_inputs import ModelSeparatedInputs

        model = ModelSeparatedInputs(
            name=run_name,
            optimizer=args.optimizer,
            learning_rate=args.learning_rate,
            regularization=args.regularization,
            dropout_conv=args.dropout_conv,
            dropout_dense=args.dropout_dense,
        )
    elif args.model == "nettcr_selu":
        from src.models.model_separated_inputs_nettcr_selu import (
            ModelSeparatedInputsNetTcrSelu,
        )

        model = ModelSeparatedInputsNetTcrSelu(
            name=run_name,
            optimizer=args.optimizer,
            learning_rate=args.learning_rate,
            regularization=args.regularization,
            dropout_conv=args.dropout_conv,
            dropout_dense=args.dropout_dense,
            activation_function=args.activation_function,
        )

    elif args.model == "nettcr":
        from src.models.model_separated_inputs_nettcr import ModelSeparatedInputsNetTcr

        model = ModelSeparatedInputsNetTcr(
            name=run_name,
            optimizer=args.optimizer,
            learning_rate=args.learning_rate,
            regularization=args.regularization,
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
            cv_type=args.cross_validation,
            n_folds=args.n_folds,
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

        train_data = separated_input_dataset_generator(
            data_stream=train,
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
            neg_shuffle=neg_shuffle,
            export_path=train_fold_output,
            neg_augment=args.neg_augment,
            augment_amount=args.augment_amount,
        )
        val_data = separated_input_dataset_generator(
            data_stream=val,
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
            neg_shuffle=neg_shuffle,
            export_path=test_fold_output,
            neg_augment=args.neg_augment,
            augment_amount=args.augment_amount,
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

        model_instance = trainer.train(model, train_data, val_data, iteration=iteration)

        trainer.clear_session(model_instance)
