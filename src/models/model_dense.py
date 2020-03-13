""" CNN model for recognizing generated peptides. """
from tensorflow.keras.layers import (
    BatchNormalization,
    Dense,
    Flatten,
)  # , Input, Dropout
from tensorflow.keras.models import Sequential

from src.models.model import Model


NUM_CLASSES = 1


class ModelDense(Model):
    def __init__(self, width, height, channels, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height
        self.channels = channels

    def _build_model(self):
        model = Sequential()

        input_shape = (self.width, self.height, self.channels)

        # model.add(Input(shape=input_shape))

        model.add(Flatten(input_shape=input_shape))
        model.add(Dense(512, activation="tanh"))
        model.add(BatchNormalization())
        model.add(Dense(512, activation="tanh"))
        model.add(BatchNormalization())
        model.add(Dense(256, activation="tanh"))
        model.add(BatchNormalization())
        model.add(Dense(NUM_CLASSES, activation="sigmoid"))

        return model

    def get_loss(self):
        from tensorflow.keras.losses import BinaryCrossentropy

        return BinaryCrossentropy()

    def get_optimizer(self):
        from tensorflow.keras.optimizers import RMSprop

        return RMSprop()
