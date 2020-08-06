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


class ModelSeparatedInputs(Model):
    def __init__(
        self,
        optimizer: str,
        learning_rate: Optional[float] = None,
        regularization: Optional[float] = None,
        dropout_conv: Optional[float] = None,
        dropout_dense: Optional[float] = None,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.optimizer = optimizer
        self.learning_rate = learning_rate
        self.regularization = l2(regularization) if regularization else None
        self.dropout_conv = dropout_conv
        self.dropout_dense = dropout_dense

    def _build_model(self):
        KERNEL_INIT = tensorflow.keras.initializers.he_normal

        input1 = Input(shape=(None, 20))
        input2 = Input(shape=(None, 20))

        def feature_extraction(input):  # noqa: A002
            convolutions = list()
            convolutions.append(
                Conv1D(
                    100,
                    (1,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="selu",
                    kernel_regularizer=self.regularization,
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (3,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="selu",
                    kernel_regularizer=self.regularization,
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (5,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="selu",
                    kernel_regularizer=self.regularization,
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (7,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="selu",
                    kernel_regularizer=self.regularization,
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (9,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="selu",
                    kernel_regularizer=self.regularization,
                )(input)
            )

            # Place them after each other (depth axis = channels = layer 2)
            merged = tensorflow.keras.layers.concatenate(convolutions, axis=-1)

            # Add dropout and batch norm (not in paper)
            if self.dropout_conv:
                merged = Dropout(self.dropout_conv)(merged)
            merged = BatchNormalization()(merged)
            conv = Conv1D(
                100,
                (1,),
                padding="valid",
                kernel_initializer=KERNEL_INIT(),
                activation="selu",
                kernel_regularizer=self.regularization,
            )(merged)
            return conv

        part1 = feature_extraction(input1)
        part2 = feature_extraction(input2)

        # axis 1 = length = sequence concat
        concatenated = tensorflow.keras.layers.concatenate([part1, part2], axis=1)

        max_pool = GlobalMaxPooling1D()(concatenated)
        if self.dropout_dense:
            drop_out = Dropout(self.dropout_dense)(max_pool)
            batch_norm = BatchNormalization()(drop_out)
        else:
            batch_norm = BatchNormalization()(max_pool)

        dense = Dense(10, activation="selu", kernel_regularizer=self.regularization,)(
            batch_norm
        )
        if self.dropout_dense:
            drop_out = Dropout(0.4)(dense)
            predictions = Dense(1, activation="sigmoid")(drop_out)
        else:
            predictions = Dense(1, activation="sigmoid")(dense)

        model = tensorflow.keras.Model(inputs=[input1, input2], outputs=predictions)

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
