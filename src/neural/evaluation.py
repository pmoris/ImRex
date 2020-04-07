import logging
import multiprocessing
from typing import Tuple

import numpy as np
import pandas as pd
import tensorflow as tf

# from tensorflow.keras import callbacks

from src.bio.feature_builder import FeatureBuilder
from src.data.vdjdb_source import VdjdbSource
from src.processing.data_stream import DataStream
from src.processing.padded_dataset_generator import padded_dataset_generator
from src.processing.separated_input_dataset_generator import (
    separated_input_dataset_generator,
)
from src.processing.zipper import Zipper


# def create_csv_logger(output_dir: Path):
#     output_path = output_dir / "metrics.csv"
#     return callbacks.CSVLogger(output_path)


# class PredictionCallBack(callbacks.Callback):
#     def __init__(self, data, base_name, output_dir: Path):
#         super().__init__()
#         self.data = data
#         self.base_name = base_name
#         self.output_dir = output_dir

#     def on_test_end(self, logs=None):
#         y_pred = self.model.predict(self.data)
#         y_true = np.array(list(self.data.unbatch().as_numpy_iterator()))[:, 1]

#         output_path = self.output_dir / "predictions.csv"

#         pd.DataFrame({"y_pred": y_pred.squeeze(), "y_true": y_true.squeeze()}).to_csv(
#             output_path, index=False
#         )


def predictions_model(model: tf.keras.Model, dataset: tf.data.Dataset):
    logger = logging.getLogger(__name__)

    logger.info("Saving predictions")
    workers = multiprocessing.cpu_count()
    logger.info(f"Using {workers} workers")

    y_pred = model.predict(dataset)
    y_true = np.array(list(dataset.unbatch().as_numpy_iterator()))[:, 1]

    predictions_df = pd.DataFrame(
        {"y_pred": y_pred.squeeze(), "y_true": y_true.squeeze()}
    )
    return predictions_df


def evaluate_model(model: tf.keras.Model, dataset: tf.data.Dataset):
    logger = logging.getLogger(__name__)

    # callbacks_list = [
    #     create_csv_logger(output_dir=output_dir),
    #     # create_tensorboard_callback(model.base_name, iteration),
    #     PredictionCallBack(
    #         data=dataset, base_name=model.base_name, output_dir=output_dir
    #     ),
    # ]

    logger.info("Evaluating CNN")
    workers = multiprocessing.cpu_count()
    logger.info(f"Using {workers} workers")

    # optionally change metrics by calling
    # model.compile(
    #     loss=model.loss, metrics=[metrics_list],
    # )

    metrics = model.evaluate(
        x=dataset,
        verbose=1,
        # callbacks=callbacks_list,
        workers=workers,
        use_multiprocessing=True,
    )

    metrics_df = pd.DataFrame([metrics], columns=model.metrics_names)
    return metrics_df


def evaluate_per_epitope(
    model: tf.keras.Model,
    data_source: VdjdbSource,
    batch_size: int,
    model_type: str,
    feature_builder: FeatureBuilder,
    cdr3_range: Tuple[int, int],
    epitope_range: Tuple[int, int],
):
    logger = logging.getLogger(__name__)

    logger.info("Evaluating CNN per epitope")
    workers = multiprocessing.cpu_count()
    logger.info(f"Using {workers} workers")

    df = pd.DataFrame()

    for epitope in data_source.data[data_source.headers["epitope_header"]].unique():

        epitope_df = data_source.data.loc[
            data_source.data[data_source.headers["epitope_header"]] == epitope
        ]

        if epitope_df["y"].unique().size != 2:
            logger.info(
                f"Skipping {epitope} because there are no negative sequence pairs."
            )
            continue

        # convert dataframe back into appropriate DataStream format for dataset generator
        test_stream = DataStream(
            Zipper(
                DataStream(zip(epitope_df["cdr3"], epitope_df["antigen.epitope"])),
                DataStream(epitope_df["y"]),
            )
        )

        if model_type == "padded":
            val_data = padded_dataset_generator(
                data_stream=test_stream,
                feature_builder=feature_builder,
                cdr3_range=cdr3_range,
                epitope_range=epitope_range,
                neg_shuffle=False,
                export_path=None,
            )
        elif model_type == "separated":
            val_data = separated_input_dataset_generator(
                data_stream=test_stream,
                cdr3_range=cdr3_range,
                epitope_range=epitope_range,
                neg_shuffle=False,
                export_path=None,
            )

        val_data = val_data.batch(batch_size)

        metrics = model.evaluate(
            x=val_data, verbose=1, workers=workers, use_multiprocessing=True,
        )
        metrics_df = pd.DataFrame([metrics], columns=model.metrics_names)
        metrics_df["epitope"] = epitope

        metrics_df["pos_data"] = epitope_df.groupby("y").size()[1]
        metrics_df["neg_data"] = epitope_df.groupby("y").size()[0]

        df = df.append(metrics_df)

    return df
