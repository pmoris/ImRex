""" Scenario for neural network with dual sequence input. """
import logging

import src.bacli as bacli
from src.config import PROJECT_ROOT
from src.data.control_cdr3_source import ControlCDR3Source
from src.data.vdjdb_source import VdjdbSource
from src.model_scripts import pipeline
from src.models.model_separated_inputs import ModelSeparatedInputs
from src.neural.trainer import get_output_path, Trainer
from src.processing.cv_folds import cv_splitter
from src.processing.separated_input_dataset_generator import (
    separated_input_dataset_generator,
)
from src.processing.splitter import splitter

bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 128,
    epochs: int = 40,
    neg_ref: bool = False,  # whether to generate negatives from CDR3 reference sequences or by shuffling positive examples
    neg_ratio: float = 0.5,  # proportion of positive to negative samples.
    val_split: float = None,  # the proportion of the dataset to include in the test split.
    epitope_grouped_cv: bool = False,  # when val_split is None, indicates whether to use normal k-fold cv or an epitope-grouped cv
    one_epitope_out_cv: bool = False,  # when val_split is None and epitope_grouped_cv is True, indicates whether to use leave-1-epitope-out cv
    neg_shuffle_in_cv: bool = True,  # NOT USED
    n_folds: int = 5,
    # these lengths are used for both size filtering and padding. Should be compatible with any preprocessing steps.
    min_length_cdr3: int = 10,
    max_length_cdr3: int = 20,
    min_length_epitope: int = 8,
    max_length_epitope: int = 13,
    name: str = "",  # name under which the model and log files will be stored, appended with the date-time.
    early_stop: bool = False,  # whether to terminate model training when the performance metric stops improving (tf.keras.callbacks.EarlyStopping)
    include_learning_rate_reduction: bool = False,  # whether to reduce the learning rate when the performance metric has stopped improving (tf.keras.callbacks.ReduceLROnPlateau)
    data_path=PROJECT_ROOT
    / "data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x.csv",  # input data csv, as supplied by preprocess_vdjdb
    optimizer: str = "rmsprop",  # can be any of: rmsprop, adam or SGD
    learning_rate: float = None,  # learning rate supplied to the selected optimizer
):

    # create logger and log file
    run_name = pipeline.create_run_name(name)
    pipeline.create_logger(run_name)
    logger = logging.getLogger(__name__)

    # log utilised function arguments that were used
    logger.info(locals())

    # read (positive) data
    data_source = VdjdbSource(
        filepath=data_path,
        headers={"cdr3_header": "cdr3", "epitope_header": "antigen.epitope"},
    )

    # check argument compatability
    if epitope_grouped_cv and val_split is not None:
        raise RuntimeError("Can't test epitope-grouped without k folds")

    logger.info("neg_ref: " + str(neg_ref))
    logger.info("epitope_grouped_cv: " + str(epitope_grouped_cv))
    logger.info("neg_shuffle_in_cv: " + str(neg_shuffle_in_cv))

    # store range restrictions for cdr3 and epitope
    cdr3_range = (min_length_cdr3, max_length_cdr3)
    epitope_range = (min_length_epitope, max_length_epitope)
    logger.info(f"cdr3 range restrictions: {cdr3_range}")
    logger.info(f"epitope range restrictions: {epitope_range}")

    trainer = Trainer(
        epochs,
        include_learning_rate_reduction=include_learning_rate_reduction,
        include_early_stop=early_stop,
    )

    model = ModelSeparatedInputs(
        name=run_name, optimizer=optimizer, learning_rate=learning_rate
    )
    logger.info(f"Built model {model.base_name}:")
    # model.summary() is logged inside trainer.py

    # if a fixed train-test split ratio is provided...
    if val_split is not None:
        train, val = splitter(data_source, test_size=val_split)
        # if a negative reference set is provided, use it
        if neg_ref:
            negative_source = ControlCDR3Source(
                min_length=min_length_cdr3, max_length=max_length_cdr3
            )
            neg_train, neg_val = splitter(negative_source, test_size=val_split)
            iterations = [((train, neg_train), (val, neg_val))]
        else:
            iterations = [(train, val)]

    # ...otherwise use a cross validation scheme
    else:
        iterations = cv_splitter(
            data_source=data_source,
            n_folds=n_folds,
            epitope_grouped=epitope_grouped_cv,
            one_out=one_epitope_out_cv,
        )

        # if a negative reference set is provided, use it
        if neg_ref:
            negative_source = ControlCDR3Source(
                min_length=min_length_cdr3, max_length=max_length_cdr3
            )
            neg_ref_fold_path = run_name + "_cdr3_ref"
            neg_iterations = cv_splitter(
                data_source=negative_source,
                n_folds=n_folds,
                epitope_grouped=False,  # CDR3 reference cannot be grouped on epitope
                one_epitope_out_cv=False,
            )
            iterations = [
                ((train, neg_train), (val, neg_val))
                for (train, val), (neg_train, neg_val) in zip(
                    iterations, neg_iterations
                )
            ]

    for iteration, (train, val) in enumerate(iterations):

        neg_train, neg_val = None, None

        logger.info(f"Iteration: {iteration}")
        logger.info(f"batch size: {batch_size}")
        if neg_ref:
            train, neg_train = train
            val, neg_val = val
            logger.info(f"neg train set: {len(neg_train)}")
            logger.info(f"neg val set: {len(neg_val)}")
        logger.info(f"train set: {len(train)}")
        logger.info(f"val set: {len(val)}")

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
            negative_ref_stream=neg_train,
            export_path=train_fold_output,
        )
        val_data = separated_input_dataset_generator(
            data_stream=val,
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
            negative_ref_stream=neg_val,
            export_path=test_fold_output,
        )

        # shuffle and batch train data
        train_data = train_data.shuffle(
            # buffer equals size of dataset, because positives and negatives are grouped
            buffer_size=train_data.reduce(0, lambda x, _: x + 1).numpy(),
            seed=42,
            # reshuffle to make each epoch see a different order of examples
            reshuffle_each_iteration=True,
        ).batch(batch_size)

        # batch validation data
        val_data = val_data.batch(batch_size)

        trainer.train(model, train_data, val_data, iteration=iteration)
