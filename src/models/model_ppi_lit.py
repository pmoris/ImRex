""" CNN model for recognizing generated peptides. """
# import keras_metrics
# from keras.layers.normalization import BatchNormalization
# from keras.regularizers import l2
import keras
import keras.initializers
from keras.layers import (
    Activation,
    Conv1D,
    Dense,
    # Dropout,
    Embedding,
    # Flatten,
    LSTM,
    MaxPool1D,
)

from src.models.model import Model


class ModelPPILit(Model):
    def __init__(self, width, height, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height

    def _build_model(self):
        # WEIGHT_DECAY = 1e-6
        KERNEL_INIT = keras.initializers.he_normal

        input1 = keras.Input(shape=(self.width,))
        input2 = keras.Input(shape=(self.height,))

        def feature_extraction(input):  # noqa: A002
            embedding = Embedding(21, 128)(input)
            x = Conv1D(10, 10, padding="valid", kernel_initializer=KERNEL_INIT())(
                embedding
            )
            x = Activation("relu")(x)
            x = MaxPool1D(pool_size=2, padding="valid")(x)
            x = Conv1D(10, 8, padding="valid", kernel_initializer=KERNEL_INIT())(x)
            x = Activation("relu")(x)
            x = MaxPool1D(pool_size=2, padding="valid")(x)
            x = Conv1D(10, 5, padding="valid", kernel_initializer=KERNEL_INIT())(x)
            x = Activation("relu")(x)
            x = MaxPool1D(pool_size=2, padding="valid")(x)
            out = LSTM(80)(x)
            # out = Flatten()(x)
            return out

        part1 = feature_extraction(input1)
        part2 = feature_extraction(input2)

        merged_vector = keras.layers.concatenate([part1, part2], axis=-1)
        predictions = Dense(1, activation="sigmoid")(merged_vector)

        model = keras.Model(inputs=[input1, input2], outputs=predictions)
        return model

    def get_loss(self):
        from keras.metrics import binary_crossentropy

        return binary_crossentropy

    def get_optimizer(self):
        from keras.optimizers import adam

        return adam()
