""" Scenario for neural network. """
import logging

import src.bacli as bacli
from src.config import PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.model_scripts import pipeline
from src.models.model_separated_inputs import ModelSeparatedInputs
from src.neural.trainer import Trainer
from src.processing.cv_folds import cv_splitter
from src.processing.separated_input_batch_generator import (
    separated_input_batch_generator,
)
from src.processing.splitter import splitter

bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 128,
    epochs: int = 40,
    neg_ratio: float = 0.5,  # proportion of positive to negative samples.
    val_split: float = None,  # the proportion of the dataset to include in the test split.
    epitope_grouped_cv: bool = False,  # when val_split is None, indicates whether to use normal k-fold cv or an epitope-grouped cv
    one_epitope_out_cv: bool = False,  # when val_split is None and epitope_grouped_cv is True, indicates whether to use leave-1-epitope-out cv
    neg_shuffle_in_cv: bool = True,  # NOT USED
    n_folds: int = 5,
    name: str = "",  # name under which the model and log files will be stored, appended with the date-time.
    min_group: int = 32,
    early_stop: bool = False,  # whether to terminate model training when the performance metric stops improving (tf.keras.callbacks.EarlyStopping)
    include_learning_rate_reduction: bool = False,  # whether to reduce the learning rate when the performance metric has stopped improving (tf.keras.callbacks.ReduceLROnPlateau)
    data_path=PROJECT_ROOT
    / "data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x.csv",  # input data csv, as supplied by preprocess_vdjdb
    optimizer: str = "rmsprop",  # can be any of: rmsprop, adam or SGD
    learning_rate: float = None,  # learning rate supplied to the selected optimizer
):

    # create logger and log file
    run_name = pipeline.create_run_name(name)
    pipeline.create_logger(name)
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

    # logger.info("neg_ref: " + str(neg_ref))
    logger.info("epitope_grouped_cv: " + str(epitope_grouped_cv))
    logger.info("neg_shuffle_in_cv: " + str(neg_shuffle_in_cv))

    trainer = Trainer(
        epochs,
        include_learning_rate_reduction=include_learning_rate_reduction,
        include_early_stop=early_stop,
    )

    model = ModelSeparatedInputs(
        optimizer=optimizer, learning_rate=learning_rate, name=name
    )
    logger.info(f"Built model {model.base_name}:")
    # model.summary() is logged inside trainer.py

    # if a fixed train-test split ratio is provided...
    if val_split is not None:
        train, val = splitter(data_source, test_size=val_split)
        iterations = [(train, val)]
    # ...otherwise use a cross validation scheme
    else:
        iterations = cv_splitter(
            data_source=data_source,
            n_folds=n_folds,
            epitope_grouped=epitope_grouped_cv,
            run_name=run_name,
        )

    for iteration, (train, val) in enumerate(iterations):
        logger.info(f"Iteration: {iteration}")
        logger.info(f"batch size: {batch_size}")
        logger.info(f"train set: {len(train)}")
        logger.info(f"val set: {len(val)}")

        train_stream = separated_input_batch_generator(
            data_stream=train,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            min_amount=min_group,
        )
        val_stream = separated_input_batch_generator(
            data_stream=val,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            min_amount=min_group,
        )

        trainer.train(model, train_stream, val_stream, iteration=iteration)
