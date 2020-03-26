""" CNN model for recognizing generated peptides. """
# import keras_metrics
# from keras.layers.normalization import BatchNormalization
# from keras.regularizers import l2
import tensorflow
from tensorflow.keras.layers import (
    Activation,
    Conv1D,
    Dense,
    Embedding,
    # Dropout,
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

        input1 = tensorflow.keras.Input(shape=(self.width,))
        input2 = tensorflow.keras.Input(shape=(self.height,))

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

        merged_vector = tensorflow.keras.layers.concatenate([part1, part2], axis=-1)
        predictions = Dense(1, activation="sigmoid")(merged_vector)

        model = tensorflow.keras.Model(inputs=[input1, input2], outputs=predictions)
        return model

    def get_loss(self):
        from tensorflow.keras.losses import BinaryCrossentropy

        return BinaryCrossentropy()

    def get_optimizer(self):
        from tensorflow.keras.optimizers import Adam

        return Adam()
