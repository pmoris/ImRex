""" CNN model for recognizing generated peptides. """
# import keras
from keras.layers import Dense, Flatten  # , Input, Dropout
from keras.layers.normalization import BatchNormalization
from keras.models import Sequential

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
        from keras.metrics import binary_crossentropy

        return binary_crossentropy

    def get_optimizer(self):
        from keras.optimizers import rmsprop

        return rmsprop()
