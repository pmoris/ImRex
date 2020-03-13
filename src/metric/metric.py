from sklearn.metrics import balanced_accuracy_score, roc_auc_score
import tensorflow as tf
from tensorflow.keras.backend import mean


def mean_pred(y_true, y_pred):
    return mean(y_pred)


def auc(y_true, y_pred):
    return tf.compat.v1.py_func(roc_auc_score, (y_true, y_pred), tf.double)


def balanced_accuracy(y_true, y_pred):
    return tf.compat.v1.py_func(
        balanced_accuracy_score, (y_true, tf.round(y_pred)), tf.double
    )
