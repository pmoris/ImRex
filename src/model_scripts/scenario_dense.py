""" Scenario for neural network. """
import logging

import src.bacli as bacli
from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.model_scripts import pipeline
from src.models.model_padded import ModelDense
from src.neural.trainer import Trainer
from src.processing.cv_folds import cv_splitter
from src.processing.inverse_map import InverseMap
from src.processing.padded_batch_generator import padded_batch_generator
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
    n_folds: int = 5,
    # these lengths are used for both size filtering and padding. Should be compatible with any preprocessing steps.
    min_length_cdr3: int = 10,
    max_length_cdr3: int = 20,
    min_length_epitope: int = 8,
    max_length_epitope: int = 13,
    name: str = "",  # name under which the model and log files will be stored, appended with the date-time.
    features: str = "hydrophob,isoelectric,mass,hydrophil,charge",  # can be any str listed in peptide_feature.featuresMap
    operator: str = "absdiff",  # can be: prod, diff, absdiff, layer or best
    early_stop: bool = False,  # whether to terminate model training when the performance metric stops improving (tf.keras.callbacks.EarlyStopping)
    include_learning_rate_reduction: bool = False,  # whether to reduce the learning rate when the performance metric has stopped improving (tf.keras.callbacks.ReduceLROnPlateau)
    data_path=PROJECT_ROOT
    / "data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x.csv",  # input data csv, as supplied by preprocess_vdjdb
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

    # get list of features and operator based on input arguments
    features_list = parse_features(features)
    operator = parse_operator(operator)
    feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)

    # check argument compatability
    if epitope_grouped_cv and val_split is not None:
        raise RuntimeError("Can't test epitope-grouped without k folds")

    logger.info("features: " + str(features_list))
    logger.info("operator: " + str(operator))
    logger.info("epitope_grouped_cv: " + str(epitope_grouped_cv))

    inverse_map = InverseMap()

    # store range restrictions for cdr3 and epitope
    cdr3_range = (min_length_cdr3, max_length_cdr3)
    epitope_range = (min_length_epitope, max_length_epitope)
    logger.info(f"cdr3 range restrictions: {cdr3_range}")
    logger.info(f"epitope range restrictions: {epitope_range}")

    trainer = Trainer(
        epochs,
        include_learning_rate_reduction=include_learning_rate_reduction,
        include_early_stop=early_stop,
        lookup=inverse_map,
    )

    model = ModelDense(
        max_length_cdr3,
        max_length_epitope,
        name=name,
        channels=feature_builder.get_number_layers(),
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
            one_out=one_epitope_out_cv,
            run_name=run_name,
        )

    for iteration, (train, val) in enumerate(iterations):
        logger.info(f"Iteration: {iteration}")
        logger.info(f"batch size: {batch_size}")
        logger.info(f"train set: {len(train)}")
        logger.info(f"val set: {len(val)}")

        train_data = padded_batch_generator(
            data_stream=train,
            feature_builder=feature_builder,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
        )
        val_data = padded_batch_generator(
            data_stream=val,
            feature_builder=feature_builder,
            neg_ratio=neg_ratio,
            batch_size=batch_size,
            cdr3_range=cdr3_range,
            epitope_range=epitope_range,
            inverse_map=inverse_map,
        )

        trainer.train(model, train_data, val_data, iteration=iteration)
