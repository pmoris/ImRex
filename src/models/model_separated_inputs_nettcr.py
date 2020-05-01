""" CNN model for recognizing generated peptides. """
from typing import Optional

# from keras import Model
# from keras.layers import LeakyReLU
# from keras.layers import Activation
# from keras.regularizers import l2
import tensorflow
from tensorflow.keras.layers import (
    BatchNormalization,
    Conv1D,
    Dense,
    Dropout,
    GlobalMaxPooling1D,
    Input,
)
from tensorflow.keras.regularizers import l2

from src.models.model import Model

NUM_CLASSES = 1
LENGTH = 10


class ModelSeparatedInputsNetTcr(Model):
    def __init__(
        self, optimizer, learning_rate: Optional[float] = None, *args, **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.optimizer = optimizer
        self.learning_rate = learning_rate

    def _build_model(self):
        n_filters = 100
        n_hid = 10
        drop_rate = 0.4

        # input:
        l_in_tcr = Input(shape=(None, 20))
        l_in_pep = Input(shape=(None, 20))

        # convolutional layers on peptide:
        l_conv_pep_1 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=1,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_pep)
        l_conv_pep_3 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=3,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_pep)
        l_conv_pep_5 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=5,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_pep)
        l_conv_pep_7 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=7,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_pep)
        l_conv_pep_9 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=9,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_pep)

        # convolutional layers on TCR:
        l_conv_tcr_1 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=1,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_tcr)
        l_conv_tcr_3 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=3,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_tcr)
        l_conv_tcr_5 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=5,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_tcr)
        l_conv_tcr_7 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=7,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_tcr)
        l_conv_tcr_9 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=9,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_in_tcr)

        # second convolutional layer:
        l_conc_1 = tensorflow.concat([l_conv_pep_1, l_conv_tcr_1], axis=1)
        l_conc_3 = tensorflow.concat([l_conv_pep_3, l_conv_tcr_3], axis=1)
        l_conc_5 = tensorflow.concat([l_conv_pep_5, l_conv_tcr_5], axis=1)
        l_conc_7 = tensorflow.concat([l_conv_pep_7, l_conv_tcr_7], axis=1)
        l_conc_9 = tensorflow.concat([l_conv_pep_9, l_conv_tcr_9], axis=1)

        l_conv_2_1 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=1,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_conc_1)
        l_conv_2_3 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=1,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_conc_3)
        l_conv_2_5 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=1,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_conc_5)
        l_conv_2_7 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=1,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_conc_7)
        l_conv_2_9 = tensorflow.keras.layers.Conv1D(
            filters=n_filters,
            kernel_size=1,
            kernel_initializer=tensorflow.keras.initializers.GlorotUniform(),
            padding="same",
            activation="sigmoid",
        )(l_conc_9)

        # max pooling:
        l_pool_max_1 = tensorflow.reduce_max(l_conv_2_1, axis=1)
        l_pool_max_3 = tensorflow.reduce_max(l_conv_2_3, axis=1)
        l_pool_max_5 = tensorflow.reduce_max(l_conv_2_5, axis=1)
        l_pool_max_7 = tensorflow.reduce_max(l_conv_2_7, axis=1)
        l_pool_max_9 = tensorflow.reduce_max(l_conv_2_9, axis=1)

        # concatenate:
        l_conc = tensorflow.concat(
            [l_pool_max_1, l_pool_max_3, l_pool_max_5, l_pool_max_7, l_pool_max_9],
            axis=1,
        )

        # dense hidden layer:
        l_dense = tensorflow.keras.layers.Dense(units=n_hid, activation="sigmoid")(
            l_conc
        )

        # dropout:
        l_dense_drop = tensorflow.keras.layers.Dropout(
            rate=drop_rate, noise_shape=None, seed=None
        )(l_dense)

        # output layer:
        out = tensorflow.keras.layers.Dense(units=1, activation="sigmoid")(l_dense_drop)

        model = tensorflow.keras.Model(inputs=[l_in_tcr, l_in_pep], outputs=out)

        return model

    def get_loss(self):
        from tensorflow.keras.losses import BinaryCrossentropy

        return BinaryCrossentropy()

    def get_optimizer(self):
        if self.optimizer == "rmsprop":

            from tensorflow.keras.optimizers import RMSprop

            if self.learning_rate:
                return RMSprop(learning_rate=self.learning_rate)
            else:
                return RMSprop()

        elif self.optimizer == "adam":

            from tensorflow.keras.optimizers import Adam

            if self.learning_rate:
                return Adam(learning_rate=self.learning_rate)
            else:
                return Adam()

        elif self.optimizer == "SGD":

            from tensorflow.keras.optimizers import SGD

            if self.learning_rate:
                return SGD(llearning_rate=self.learning_rate)
            else:
                return SGD()
