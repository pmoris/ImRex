import keras.backend as K
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
import tensorflow as tf


def mean_pred(y_true, y_pred):
    return K.mean(y_pred)


def AUC(y_true, y_pred):
    return tf.py_func(roc_auc_score, (y_true, y_pred), tf.double)


def balanced_accuracy(y_true, y_pred):
    return tf.py_func(balanced_accuracy_score, (y_true, tf.round(y_pred)), tf.double)
