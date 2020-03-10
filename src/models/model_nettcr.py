""" CNN model for recognizing generated peptides. """
# from keras import Model
# from keras.layers import LeakyReLU
# from keras.layers import Activation
# from keras.regularizers import l2
import keras
import keras.initializers
from keras.layers import (
    Conv1D,
    Dense,
    # Dropout,
    # Flatten,
    # Conv2D,
    # MaxPool2D,
    Input,
)
from keras.layers import (
    # GlobalAveragePooling2D,
    # GlobalMaxPooling2D,
    GlobalMaxPooling1D,
)
from keras.layers.normalization import BatchNormalization

from src.models.model import Model

NUM_CLASSES = 1
LENGTH = 10


class ModelNetTCR(Model):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

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
                    activation="sigmoid",
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (3,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="sigmoid",
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (5,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="sigmoid",
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (7,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="sigmoid",
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (9,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT(),
                    activation="sigmoid",
                )(input)
            )

            # Place them after each other (depth axis = channels = layer 2)
            merged = keras.layers.concatenate(convolutions, axis=-1)

            # Add dropout and batch norm (not in paper)
            # merged = Dropout(0.25)(merged)
            merged = BatchNormalization()(merged)
            conv = Conv1D(
                100,
                (1,),
                padding="valid",
                kernel_initializer=KERNEL_INIT(),
                activation="sigmoid",
            )(merged)
            return conv

        part1 = feature_extraction(input1)
        part2 = feature_extraction(input2)

        concatenated = keras.layers.concatenate([part1, part2], axis=1)

        max_pool = GlobalMaxPooling1D()(concatenated)

        dense = Dense(10, activation="sigmoid")(max_pool)
        # dense = Dropout(0.5)(dense)
        predictions = Dense(1, activation="sigmoid")(dense)

        model = keras.Model(inputs=[input1, input2], outputs=predictions)
        return model

    def get_loss(self):
        from keras.metrics import binary_crossentropy

        return binary_crossentropy

    def get_optimizer(self):
        from keras.optimizers import adam

        return adam()
