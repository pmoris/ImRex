import os
from pathlib import Path

from keras.callbacks import Callback
from keras.utils import multi_gpu_model
import keras_metrics
import multiprocessing
import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_curve,
    auc,
    precision_recall_curve,
    average_precision_score,
)

from src.metric import metric
from src.config import PROJECT_ROOT


NUMBER_OF_GPUS = int(os.environ["GPUS"])

LOGDIR = PROJECT_ROOT / "models/logs"
OUTDIR = PROJECT_ROOT / "models/models"


def get_output_dir(base_name, iteration=None):
    if iteration is None:
        ret = os.path.join(OUTDIR, base_name)
    else:
        ret = os.path.join(OUTDIR, base_name, "iteration {}".format(iteration))
    os.makedirs(ret, exist_ok=True)
    return ret


def get_output_path(base_name, file_name, iteration=None):
    return os.path.join(get_output_dir(base_name, iteration), file_name)


def create_checkpointer(base_name, iteration):
    from keras.callbacks import ModelCheckpoint

    output_path = get_output_path(base_name, base_name + ".h5", iteration=iteration)
    return ModelCheckpoint(
        output_path, monitor="val_loss", verbose=1, save_best_only=True, mode="min"
    )


def create_tensorboard_callback(model_name):
    from keras.callbacks import TensorBoard

    logpath = os.path.join(LOGDIR, model_name)

    return TensorBoard(log_dir=logpath)


def create_early_stopping(patience=5):
    from keras.callbacks import EarlyStopping

    return EarlyStopping(patience=patience)


def create_csv_logger(base_name, iteration):
    from keras.callbacks import CSVLogger

    output_path = get_output_path(base_name, "metrics.csv", iteration=iteration)
    return CSVLogger(output_path)


def create_LRR():
    from keras.callbacks import ReduceLROnPlateau

    return ReduceLROnPlateau(
        monitor="loss", patience=3, verbose=1, factor=0.5, min_lr=0.0001
    )


class MetricCallback(Callback):
    def __init__(self, val_stream, base_name, iteration):
        super().__init__()
        self.val_stream = val_stream
        self.base_name = base_name
        self.iteration = iteration

    def _predict(self, include_samples=False):
        xx = list()
        y_pred = list()
        y_true = list()
        for i in range(len(self.val_stream) * 5):
            batch = self.val_stream[i]
            x, y = batch
            pred = self.model.predict_on_batch(x)
            if include_samples:
                xx.extend(x)
            y_pred.extend(pred)
            y_true.extend(y)

        if include_samples:
            return y_pred, y_true, xx
        else:
            return y_pred, y_true


class RocCallback(MetricCallback):
    def __init__(self, val_stream, base_name, iteration):
        super().__init__(val_stream, base_name, iteration)

    def on_train_end(self, logs={}):
        y_pred, y_true = self._predict()

        fpr, tpr, _ = roc_curve(y_true, y_pred, drop_intermediate=True)
        auc_value = auc(fpr, tpr)

        interval = np.linspace(0, 1, 201)
        tpr_i = np.interp(interval, fpr, tpr)

        output_path = get_output_path(
            self.base_name, "roc.csv", iteration=self.iteration
        )
        df = pd.DataFrame({"fpr": interval, "tpr": tpr_i})
        df.to_csv(output_path, index=False)

        output_path = get_output_path(
            self.base_name, "auc.csv", iteration=self.iteration
        )
        df = pd.DataFrame({"auc": [auc_value],})
        df.to_csv(output_path, index=False)

        return


class PrecisionRecallCallback(MetricCallback):
    def __init__(self, val_stream, base_name, iteration):
        super().__init__(val_stream, base_name, iteration)

    def on_train_end(self, logs={}):
        y_pred, y_true = self._predict()
        precision, recall, _ = precision_recall_curve(y_true, y_pred)
        ap = average_precision_score(y_true, y_pred)

        interval = np.linspace(0, 1, 201)
        precision_i = np.interp(interval, recall[::-1], precision[::-1])

        output_path = get_output_path(
            self.base_name, "precision_recall.csv", iteration=self.iteration
        )
        df = pd.DataFrame({"recall": interval, "precision": precision_i})
        df.to_csv(output_path, index=False)

        output_path = get_output_path(
            self.base_name, "average_precision.csv", iteration=self.iteration
        )
        df = pd.DataFrame({"average_precision": [ap],})
        df.to_csv(output_path, index=False)

        return


class PredictionCallback(MetricCallback):
    def __init__(self, val_stream, base_name, iteration, lookup=None):
        super().__init__(val_stream, base_name, iteration)
        self.lookup = lookup

    def on_train_end(self, logs={}):
        if self.lookup:
            y_pred, y_true, samples = self._predict(include_samples=True)
        else:
            y_pred, y_true = self._predict(include_samples=False)

        y_pred = [e[0] for e in y_pred]
        df = pd.DataFrame(zip(y_true, y_pred), columns=["y_true", "y_pred"])
        output_path = get_output_path(
            self.base_name, "predictions.csv", iteration=self.iteration
        )
        df.to_csv(output_path, index=False)

        if self.lookup:
            all_pos = [
                (pred, sample, 1)
                for pred, sample, label in zip(y_pred, samples, y_true)
                if label == 1
            ]
            all_neg = [
                (pred, sample, 0)
                for pred, sample, label in zip(y_pred, samples, y_true)
                if label == 0
            ]

            interesting = {
                "Best Positive": max(all_pos, key=lambda x: x[0]),
                "Worst Positive": min(all_pos, key=lambda x: x[0]),
                "Best Negative": min(all_neg, key=lambda x: x[0]),
                "Worst Negative": max(all_neg, key=lambda x: x[0]),
            }

            for name, (prediction, sample, label) in interesting.items():
                print(f" ======= {name} ======= ")
                print(f"prediction:\t{prediction}")
                print(f"label:\t\t{label}")
                x = self.lookup.find_input_for(sample)
                if x:
                    cdr3, epitope = x
                    print(f"epitope:\t{epitope}")
                    print(f"cdr3:\t\t{cdr3}")
                else:
                    print("Not found")


def get_metrics():
    return [
        "accuracy",
        metric.balanced_accuracy,
        metric.mean_pred,
        metric.AUC,
        keras_metrics.precision(),
        keras_metrics.recall(),
    ]


class Trainer(object):
    def __init__(
        self, epochs, include_LRR=False, lookup=None, include_early_stop=False
    ):
        self.epochs = epochs
        self.include_LRR = include_LRR
        self.include_early_stop = include_early_stop
        self.histories = dict()
        self.base_name = None
        self.lookup = lookup  # Lookup map from feature to input (for traceability)

    def train(self, model, train_stream, val_stream, iteration=None):
        model_instance = model.new_instance()

        # Print summary once, and before converting to multi GPU
        if iteration == 0:
            print("Training model:")
            model_instance.summary()

        if NUMBER_OF_GPUS > 1:
            model_instance = multi_gpu_model(model_instance, gpus=NUMBER_OF_GPUS)

        model_instance.compile(
            loss=model.get_loss(),
            optimizer=model.get_optimizer(),
            metrics=get_metrics(),
        )

        if not self.base_name:
            self.base_name = model.base_name

        callbacks = [
            create_checkpointer(model.base_name, iteration),
            create_csv_logger(model.base_name, iteration),
            create_tensorboard_callback(model.get_name(iteration)),
            RocCallback(val_stream, model.base_name, iteration),
            PrecisionRecallCallback(val_stream, model.base_name, iteration),
            PredictionCallback(
                val_stream, model.base_name, iteration, lookup=self.lookup
            ),
        ]

        if self.include_early_stop:
            callbacks.append(create_early_stopping())

        if self.include_LRR:
            callbacks.append(create_LRR())

        print("Fitting CNN")
        workers = multiprocessing.cpu_count()
        print(f"Using {workers} workers")
        history = model_instance.fit_generator(
            generator=train_stream,
            epochs=self.epochs,
            verbose=1,
            callbacks=callbacks,
            validation_data=val_stream,
            class_weight=None,
            max_queue_size=2,
            workers=workers,
            use_multiprocessing=True,
            shuffle=False,
        )
        self.histories[iteration] = history

    def evaluate(self):
        pass  # Moved to post processing step
