import logging
import multiprocessing
import os
from typing import Optional

import numpy as np
import pandas as pd
from tensorflow.keras import callbacks
from tensorflow.keras import metrics
from tensorflow.keras.utils import multi_gpu_model

from src.config import MODEL_DIR, TENSORBOARD_DIR
from src.processing.inverse_map import InverseMap


NUMBER_OF_GPUS = os.environ.get("GPUS")


def get_output_dir(base_name, iteration=None):
    """Return and create the output directories that will hold the model (iteration) output."""
    if iteration is None:
        ret = MODEL_DIR / base_name
    else:
        ret = MODEL_DIR / base_name / f"iteration_{iteration}"
    ret.mkdir(parents=True, exist_ok=True)
    return ret


def get_output_path(base_name, file_name, iteration=None):
    """Return the output filepath to store the model (iteration)."""
    output_path = get_output_dir(base_name, iteration) / file_name
    return str(output_path.resolve())


def create_checkpointer(base_name, iteration):
    output_path = get_output_path(
        base_name=base_name,
        file_name=base_name + "-epoch{epoch:02d}-valacc{val_accuracy:.2f}.h5",
        iteration=iteration,
    )

    return callbacks.ModelCheckpoint(
        filepath=output_path,
        verbose=1,
        # Path where to save the model
        # The two parameters below mean that we will overwrite
        # the current checkpoint if and only if
        # the `val_loss` score has improved.
        save_best_only=True,
        monitor="val_loss",
        mode="auto",  # For val_acc, this should be max, for val_loss this should be min, etc. In auto mode, the direction is automatically inferred from the name of the monitored quantity.
    )


def create_tensorboard_callback(base_name, iteration):
    logpath = TENSORBOARD_DIR / base_name / f"iteration_{iteration}"
    logpath.parent.mkdir(parents=True, exist_ok=True)

    return callbacks.TensorBoard(log_dir=logpath)


def create_csv_logger(base_name, iteration):
    output_path = get_output_path(base_name, "metrics.csv", iteration=iteration)

    return callbacks.CSVLogger(output_path)


# class MetricCallback(callbacks.Callback):
#     """Callback template to compute average metrics after training."""

#     def __init__(self, val_stream, base_name, iteration):
#         super().__init__()
#         self.val_stream = val_stream
#         self.base_name = base_name
#         self.iteration = iteration

#     def _predict(self, include_samples=False):
#         xx = list()
#         y_pred = list()
#         y_true = list()
#         for i in range(len(self.val_stream) * 5):
#             batch = self.val_stream[i]
#             x, y = batch
#             pred = self.model.predict_on_batch(x)
#             if include_samples:
#                 xx.extend(x)
#             y_pred.extend(pred)
#             y_true.extend(y)

#         if include_samples:
#             return y_pred, y_true, xx
#         else:
#             return y_pred, y_true


# class RocCallback(MetricCallback):
#     def __init__(self, val_stream, base_name, iteration):
#         super().__init__(val_stream, base_name, iteration)
#
#     def on_train_end(self, logs={}):
#         y_pred, y_true = self._predict()
#
#         fpr, tpr, _ = roc_curve(y_true, y_pred, drop_intermediate=True)
#         auc_value = auc(fpr, tpr)
#
#         interval = np.linspace(0, 1, 201)
#         tpr_i = np.interp(interval, fpr, tpr)
#
#         output_path = get_output_path(
#             self.base_name, "roc.csv", iteration=self.iteration
#         )
#         df = pd.DataFrame({"fpr": interval, "tpr": tpr_i})
#         df.to_csv(output_path, index=False)
#
#         output_path = get_output_path(
#             self.base_name, "auc.csv", iteration=self.iteration
#         )
#         df = pd.DataFrame({"auc": [auc_value]})
#         df.to_csv(output_path, index=False)
#
#         return
#
#
# class PrecisionRecallCallback(MetricCallback):
#     def __init__(self, val_stream, base_name, iteration):
#         super().__init__(val_stream, base_name, iteration)
#
#     def on_train_end(self, logs={}):
#         y_pred, y_true = self._predict()
#         precision, recall, _ = precision_recall_curve(y_true, y_pred)
#         ap = average_precision_score(y_true, y_pred)
#
#         interval = np.linspace(0, 1, 201)
#         precision_i = np.interp(interval, recall[::-1], precision[::-1])
#
#         output_path = get_output_path(
#             self.base_name, "precision_recall.csv", iteration=self.iteration
#         )
#         df = pd.DataFrame({"recall": interval, "precision": precision_i})
#         df.to_csv(output_path, index=False)
#
#         output_path = get_output_path(
#             self.base_name, "average_precision.csv", iteration=self.iteration
#         )
#         df = pd.DataFrame({"average_precision": [ap]})
#         df.to_csv(output_path, index=False)
#
#         return


# class PredictionCallback(MetricCallback):
#     def __init__(self, val_stream, base_name, iteration, lookup=None):
#         super().__init__(val_stream, base_name, iteration)
#         self.lookup = lookup

#     def on_train_end(self, logs={}):
#         if self.lookup:
#             y_pred, y_true, samples = self._predict(include_samples=True)
#         else:
#             y_pred, y_true = self._predict(include_samples=False)

#         y_pred = [e[0] for e in y_pred]
#         df = pd.DataFrame(zip(y_true, y_pred), columns=["y_true", "y_pred"])
#         output_path = get_output_path(
#             self.base_name, "predictions.csv", iteration=self.iteration
#         )
#         df.to_csv(output_path, index=False)

#         if self.lookup:
#             all_pos = [
#                 (pred, sample, 1)
#                 for pred, sample, label in zip(y_pred, samples, y_true)
#                 if label == 1
#             ]
#             all_neg = [
#                 (pred, sample, 0)
#                 for pred, sample, label in zip(y_pred, samples, y_true)
#                 if label == 0
#             ]

#             interesting = {
#                 "Best Positive": max(all_pos, key=lambda x: x[0]),
#                 "Worst Positive": min(all_pos, key=lambda x: x[0]),
#                 "Best Negative": min(all_neg, key=lambda x: x[0]),
#                 "Worst Negative": max(all_neg, key=lambda x: x[0]),
#             }

#             for name, (prediction, sample, label) in interesting.items():
#                 print(f" ======= {name} ======= ")
#                 print(f"prediction:\t{prediction}")
#                 print(f"label:\t\t{label}")
#                 x = self.lookup.find_input_for(sample)
#                 if x:
#                     cdr3, epitope = x
#                     print(f"epitope:\t{epitope}")
#                     print(f"cdr3:\t\t{cdr3}")
#                 else:
#                     print("Not found")


class PredictionCallBack(callbacks.Callback):
    def __init__(self, val_data, base_name, iteration, lookup=None):
        super().__init__()
        self.val_data = val_data
        self.base_name = base_name
        self.iteration = iteration
        self.lookup = lookup

    def on_train_end(self, logs=None):
        y_pred = self.model.predict(self.val_data)
        y_true = np.array(list(self.val_data.unbatch().as_numpy_iterator()))[:, 1]

        output_path = get_output_path(
            self.base_name, f"predictions.csv", iteration=self.iteration
        )

        # if self.lookup:
        #     samples = [
        #         self.lookup.find_input_for(sample)
        #         for sample in self.val_data.unbatch().as_numpy_iterator()
        #     ]
        #     cdr3, epitope = zip(*samples)
        #     pd.DataFrame(
        #         {
        #             "y_pred": y_pred.squeeze(),
        #             "y_true": y_true.squeeze(),
        #             "cdr3": cdr3,
        #             "antigen.epitope": epitope,
        #         }
        #     ).to_csv(output_path, index=False)
        # else:
        pd.DataFrame({"y_pred": y_pred.squeeze(), "y_true": y_true.squeeze()}).to_csv(
            output_path, index=False
        )


def get_metrics():
    return [
        metrics.BinaryAccuracy(name="accuracy"),
        metrics.AUC(curve="ROC", name="roc_auc"),
        metrics.AUC(curve="PR", name="pr_auc"),
        metrics.Precision(name="precision"),
        metrics.Recall(name="recall"),
        metrics.TruePositives(name="tp"),
        metrics.FalsePositives(name="fp"),
        metrics.TrueNegatives(name="tn"),
        metrics.FalseNegatives(name="fn"),
    ]


class Trainer(object):
    """Object that can compile and train a keras model instance with callbacks."""

    def __init__(
        self,
        epochs: int,
        include_learning_rate_reduction: bool = False,
        include_early_stop: bool = False,
        lookup: Optional[InverseMap] = None,
        verbose: bool = True,
    ):
        """Initialise Trainer object.

        Parameters
        ----------
        epochs : int
            Number of epochs to train the model, passed to keras' Model.fit().
        include_learning_rate_reduction : bool, optional
            Add callback to reduce learning rate when a metric has stopped improving, by default False
        lookup : Optional[InverseMap], optional
            [description], by default None
        include_early_stop : bool, optional
            [description], by default False
        verbose : bool, optional
            True for verbose training output, False otherwise, by default True.
        """
        self.epochs = epochs
        self.include_learning_rate_reduction = include_learning_rate_reduction
        self.include_early_stop = include_early_stop
        self.histories = dict()
        self.base_name = None
        self.lookup = lookup  # Lookup map from feature to input (for traceability)
        self.verbose = 1 if verbose else 0

    def train(self, model, train_data, val_data, iteration=None):
        logger = logging.getLogger(__name__)

        model_instance = model.new_instance()

        # Print summary once, and before converting to multi GPU
        if iteration == 0:
            logger.info("Training model:")
            model_instance.summary(print_fn=logger.info)

        if NUMBER_OF_GPUS and int(NUMBER_OF_GPUS) > 1:
            model_instance = multi_gpu_model(model_instance, gpus=int(NUMBER_OF_GPUS))

        model_instance.compile(
            optimizer=model.get_optimizer(),
            loss=model.get_loss(),
            metrics=get_metrics(),
            loss_weights=None,
            sample_weight_mode=None,
            weighted_metrics=None,
            target_tensors=None,
            distribute=None,
        )

        if not self.base_name:
            self.base_name = model.base_name

        callbacks_list = [
            create_checkpointer(model.base_name, iteration),
            create_csv_logger(model.base_name, iteration),
            create_tensorboard_callback(model.base_name, iteration),
            PredictionCallBack(val_data, model.base_name, iteration, self.lookup)
            # RocCallback(val_stream, model.base_name, iteration),
            # PrecisionRecallCallback(val_stream, model.base_name, iteration),
            # PredictionCallback(
            #     val_stream, model.base_name, iteration, lookup=self.lookup
            # ),
        ]

        if self.include_early_stop:
            callbacks_list.append(
                callbacks.EarlyStopping(  # Stop training when `val_loss` is no longer improving
                    monitor="val_loss",
                    # "no longer improving" being defined as "no better than 1e-2 less"
                    min_delta=1e-2,
                    # "no longer improving" being further defined as "for at least 2 epochs"
                    patience=10,
                    verbose=1,
                )
            )

        if self.include_learning_rate_reduction:
            callbacks_list.append(
                callbacks.ReduceLROnPlateau(
                    monitor="val_loss",
                    factor=0.2,  # factor by which the learning rate will be reduced. new_lr = lr * factor
                    patience=3,  # number of epochs with no improvement after which learning rate will be reduced
                    verbose=1,
                    mode="auto",
                    min_delta=0.0001,  # threshold for measuring the new optimum, to only focus on significant changes
                    cooldown=0,  # number of epochs to wait before resuming normal operation after lr has been reduced
                    min_lr=0,  # lower bound on the learning rate
                )
            )

        # save model with weight initialization before training
        initial_model_path = get_output_path(
            base_name=model.base_name,
            file_name=model.base_name + "-init.h5",
            iteration=iteration,
        )
        # model_instance.evaluate(val_data) # get metrics and print them to file + store for modelfilename
        model_instance.save(initial_model_path)

        logger.info("Fitting CNN")
        workers = multiprocessing.cpu_count()
        logger.info(f"Using {workers} workers")
        history = model_instance.fit(
            x=train_data,
            epochs=self.epochs,
            verbose=self.verbose,
            callbacks=callbacks_list,
            validation_data=val_data,
            class_weight=None,
            max_queue_size=2,
            workers=workers,
            use_multiprocessing=True,
            shuffle=False,
        )
        self.histories[iteration] = history

        return model_instance

    def evaluate(self):
        pass  # Moved to post processing step
