""" Scenario for neural network. """
import datetime
import logging

import src.bacli as bacli
from src.config import LOG_DIR, PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.models.model_separated_inputs import ModelSeparatedInputs
from src.neural.trainer import Trainer
from src.processing.cv_folds import cv_splitter
from src.processing.net_tcr_batch_generator import nettcr_batch_generator
from src.processing.splitter import splitter

bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 128,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    val_split: float = None,  # the proportion of the dataset to include in the test split.
    epitope_grouped_cv: bool = False,
    neg_shuffle_in_cv: bool = True,
    n_folds: int = 5,
    name: str = "",
    min_group: int = 32,
    early_stop=False,
    include_learning_rate_reduction: bool = False,
    data_path=PROJECT_ROOT
    / "data/interim/vdjdb-2020-01-20/vdjdb-human-tra-trb-no10x.csv",
    optimizer: str = "rmsprop",  # can be any of: rmsprop, adam or SGD
    learning_rate: bool = False,
):

    # create run name by appending time and date
    run_name = name + datetime.datetime.now().strftime("_%Y%m%d_%H-%M-%S")
    # create filepath for log
    log_file = LOG_DIR / run_name
    log_file = log_file.with_suffix(".log")
    log_file.parent.mkdir(parents=True, exist_ok=True)
    # create file logger
    log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(filename=log_file, level=logging.INFO, format=log_fmt)
    # apply settings to root logger, so that loggers in modules can inherit both the file and console logger
    logger = logging.getLogger()
    # add console logger
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter(log_fmt))
    logger.addHandler(console)

    # log utilised function arguments that were used for logging purposes
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
        optimizer=optimizer, learning_rate=learning_rate, name_suffix=name
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

        train_stream = nettcr_batch_generator(
            data_stream=train,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            min_amount=min_group,
        )
        val_stream = nettcr_batch_generator(
            data_stream=val,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            min_amount=min_group,
        )

        trainer.train(model, train_stream, val_stream, iteration=iteration)
