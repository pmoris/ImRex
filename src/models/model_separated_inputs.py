""" CNN model for recognizing generated peptides. """
from typing import Optional

# from keras import Model
# from keras.layers import LeakyReLU
# from keras.layers import Activation
# from keras.regularizers import l2
import keras
import keras.initializers
from keras.layers import (
    Conv1D,
    Dense,
    Dropout,
    GlobalMaxPooling1D,
    Input,
)
from keras.layers.normalization import BatchNormalization

from src.models.model import Model

NUM_CLASSES = 1
LENGTH = 10


class ModelSeparatedInputs(Model):
    def __init__(
        self, optimizer, learning_rate: Optional[bool] = None, *args, **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.optimizer = optimizer
        self.learning_rate = learning_rate

    def _build_model(self):
        KERNEL_INIT = keras.initializers.he_normal

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
                    activation="relu",
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (3,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="relu",
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (5,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="relu",
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (7,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="relu",
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (9,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="relu",
                )(input)
            )

            # Place them after each other (depth axis = channels = layer 2)
            merged = keras.layers.concatenate(convolutions, axis=-1)

            # Add dropout and batch norm (not in paper)
            merged = Dropout(0.25)(merged)
            merged = BatchNormalization()(merged)
            conv = Conv1D(
                100,
                (1,),
                padding="valid",
                kernel_initializer=KERNEL_INIT(),
                activation="relu",
            )(merged)
            return conv

        part1 = feature_extraction(input1)
        part2 = feature_extraction(input2)

        # axis 1 = length = sequence concat
        concatenated = keras.layers.concatenate([part1, part2], axis=1)

        max_pool = GlobalMaxPooling1D()(concatenated)
        drop_out = Dropout(0.5)(max_pool)
        batch_norm = BatchNormalization()(drop_out)

        dense = Dense(10, activation="relu")(batch_norm)
        predictions = Dense(1, activation="sigmoid")(dense)

        model = keras.Model(inputs=[input1, input2], outputs=predictions)

        return model

    def get_loss(self):
        from keras.metrics import binary_crossentropy

        return binary_crossentropy

    def get_optimizer(self):
        if self.optimizer == "rmsprop":

            from keras.optimizers import RMSprop

            if self.learning_rate:
                return RMSprop(learning_rate=self.learning_rate)
            else:
                return RMSprop()

        elif self.optimizer == "adam":

            from keras.optimizers import Adam

            if self.learning_rate:
                return Adam(learning_rate=self.learning_rate)
            else:
                return Adam()

        elif self.optimizer == "SGD":

            from keras.optimizers import SGD

            if self.learning_rate:
                return SGD(llearning_rate=self.learning_rate)
            else:
                return SGD()
