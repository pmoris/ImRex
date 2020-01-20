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


def getOutputDir(baseName, iteration=None):
    if iteration is None:
        ret = os.path.join(OUTDIR, baseName)
    else:
        ret = os.path.join(OUTDIR, baseName, "iteration {}".format(iteration))
    os.makedirs(ret, exist_ok=True)
    return ret


def getOutputPath(baseName, fileName, iteration=None):
    return os.path.join(getOutputDir(baseName, iteration), fileName)


def createCheckpointer(baseName, iteration):
    from keras.callbacks import ModelCheckpoint

    outputPath = getOutputPath(baseName, baseName + ".h5", iteration=iteration)
    return ModelCheckpoint(
        outputPath, monitor="val_loss", verbose=1, save_best_only=True, mode="min"
    )


def createTensorboardCallback(modelName):
    from keras.callbacks import TensorBoard

    logpath = os.path.join(LOGDIR, modelName)

    return TensorBoard(log_dir=logpath)


def createEarlyStopping(patience=5):
    from keras.callbacks import EarlyStopping

    return EarlyStopping(patience=patience)


def createCSVLogger(baseName, iteration):
    from keras.callbacks import CSVLogger

    outputPath = getOutputPath(baseName, "metrics.csv", iteration=iteration)
    return CSVLogger(outputPath)


def createLRR():
    from keras.callbacks import ReduceLROnPlateau

    return ReduceLROnPlateau(
        monitor="loss", patience=3, verbose=1, factor=0.5, min_lr=0.0001
    )


class MetricCallback(Callback):
    def __init__(self, valStream, baseName, iteration):
        super().__init__()
        self.valStream = valStream
        self.baseName = baseName
        self.iteration = iteration

    def _predict(self, includeSamples=False):
        xx = list()
        y_pred = list()
        y_true = list()
        for i in range(len(self.valStream) * 5):
            batch = self.valStream[i]
            x, y = batch
            pred = self.model.predict_on_batch(x)
            if includeSamples:
                xx.extend(x)
            y_pred.extend(pred)
            y_true.extend(y)

        if includeSamples:
            return y_pred, y_true, xx
        else:
            return y_pred, y_true


class roc_callback(MetricCallback):
    def __init__(self, valStream, baseName, iteration):
        super().__init__(valStream, baseName, iteration)

    def on_train_end(self, logs={}):
        y_pred, y_true = self._predict()

        fpr, tpr, _ = roc_curve(y_true, y_pred, drop_intermediate=True)
        aucValue = auc(fpr, tpr)

        interval = np.linspace(0, 1, 201)
        tprI = np.interp(interval, fpr, tpr)

        outputPath = getOutputPath(self.baseName, "roc.csv", iteration=self.iteration)
        df = pd.DataFrame({"fpr": interval, "tpr": tprI})
        df.to_csv(outputPath, index=False)

        outputPath = getOutputPath(self.baseName, "auc.csv", iteration=self.iteration)
        df = pd.DataFrame({"auc": [aucValue],})
        df.to_csv(outputPath, index=False)

        return


class precision_recall_callback(MetricCallback):
    def __init__(self, valStream, baseName, iteration):
        super().__init__(valStream, baseName, iteration)

    def on_train_end(self, logs={}):
        y_pred, y_true = self._predict()
        precision, recall, _ = precision_recall_curve(y_true, y_pred)
        ap = average_precision_score(y_true, y_pred)

        interval = np.linspace(0, 1, 201)
        precisionI = np.interp(interval, recall[::-1], precision[::-1])

        outputPath = getOutputPath(
            self.baseName, "precision_recall.csv", iteration=self.iteration
        )
        df = pd.DataFrame({"recall": interval, "precision": precisionI})
        df.to_csv(outputPath, index=False)

        outputPath = getOutputPath(
            self.baseName, "average_precision.csv", iteration=self.iteration
        )
        df = pd.DataFrame({"average_precision": [ap],})
        df.to_csv(outputPath, index=False)

        return


class prediction_callback(MetricCallback):
    def __init__(self, valStream, baseName, iteration, lookup=None):
        super().__init__(valStream, baseName, iteration)
        self.lookup = lookup

    def on_train_end(self, logs={}):
        if self.lookup:
            y_pred, y_true, samples = self._predict(includeSamples=True)
        else:
            y_pred, y_true = self._predict(includeSamples=False)

        y_pred = [e[0] for e in y_pred]
        df = pd.DataFrame(zip(y_true, y_pred), columns=["y_true", "y_pred"])
        outputPath = getOutputPath(
            self.baseName, "predictions.csv", iteration=self.iteration
        )
        df.to_csv(outputPath, index=False)

        if self.lookup:
            allPos = [
                (pred, sample, 1)
                for pred, sample, label in zip(y_pred, samples, y_true)
                if label == 1
            ]
            allNeg = [
                (pred, sample, 0)
                for pred, sample, label in zip(y_pred, samples, y_true)
                if label == 0
            ]

            interesting = {
                "Best Positive": max(allPos, key=lambda x: x[0]),
                "Worst Positive": min(allPos, key=lambda x: x[0]),
                "Best Negative": min(allNeg, key=lambda x: x[0]),
                "Worst Negative": max(allNeg, key=lambda x: x[0]),
            }

            for name, (prediction, sample, label) in interesting.items():
                print(f" ======= {name} ======= ")
                print(f"prediction:\t{prediction}")
                print(f"label:\t\t{label}")
                x = self.lookup.findInputFor(sample)
                if x:
                    cdr3, epitope = x
                    print(f"epitope:\t{epitope}")
                    print(f"cdr3:\t\t{cdr3}")
                else:
                    print("Not found")


def getMetrics():
    return [
        "accuracy",
        metric.balanced_accuracy,
        metric.mean_pred,
        metric.AUC,
        keras_metrics.precision(),
        keras_metrics.recall(),
    ]


class Trainer(object):
    def __init__(self, epochs, includeLRR=False, lookup=None, includeEarlyStop=False):
        self.epochs = epochs
        self.includeLRR = includeLRR
        self.includeEarlyStop = includeEarlyStop
        self.histories = dict()
        self.baseName = None
        self.lookup = lookup  # Lookup map from feature to input (for traceability)

    def train(self, model, trainStream, valStream, iteration=None):
        modelInstance = model.newInstance()

        # Print summary once, and before converting to multi GPU
        if iteration == 0:
            print("Training model:")
            modelInstance.summary()

        if NUMBER_OF_GPUS > 1:
            modelInstance = multi_gpu_model(modelInstance, gpus=NUMBER_OF_GPUS)

        modelInstance.compile(
            loss=model.getLoss(), optimizer=model.getOptimizer(), metrics=getMetrics()
        )

        if not self.baseName:
            self.baseName = model.baseName

        callbacks = [
            createCheckpointer(model.baseName, iteration),
            createCSVLogger(model.baseName, iteration),
            createTensorboardCallback(model.getName(iteration)),
            roc_callback(valStream, model.baseName, iteration),
            precision_recall_callback(valStream, model.baseName, iteration),
            prediction_callback(
                valStream, model.baseName, iteration, lookup=self.lookup
            ),
        ]

        if self.includeEarlyStop:
            callbacks.append(createEarlyStopping())

        if self.includeLRR:
            callbacks.append(createLRR())

        print("Fitting CNN")
        workers = multiprocessing.cpu_count()
        print(f"Using {workers} workers")
        history = modelInstance.fit_generator(
            generator=trainStream,
            epochs=self.epochs,
            verbose=1,
            callbacks=callbacks,
            validation_data=valStream,
            class_weight=None,
            max_queue_size=2,
            workers=workers,
            use_multiprocessing=True,
            shuffle=False,
        )
        self.histories[iteration] = history

    def evaluate(self):
        pass  # Moved to post processing step
