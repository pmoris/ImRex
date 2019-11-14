""" CNN model for recognizing generated peptides. """
import os

from .model import Model


NUM_CLASSES = 1


class ModelDense(Model):
    def __init__(self, width, height, channels, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height
        self.channels = channels

    def _buildModel(self):
        import keras
        from keras.models import Sequential
        from keras.layers import Dense, Flatten, Input, Dropout
        from keras.layers.normalization import BatchNormalization

        model = Sequential()

        inputShape = (self.width, self.height, self.channels)

        # model.add(Input(shape=inputShape))

        model.add(Flatten(input_shape=inputShape))
        model.add(Dense(512, activation='tanh'))
        model.add(BatchNormalization())
        model.add(Dense(512, activation='tanh'))
        model.add(BatchNormalization())
        model.add(Dense(256, activation='tanh'))
        model.add(BatchNormalization())
        model.add(Dense(NUM_CLASSES, activation='sigmoid'))

        return model

    def getLoss(self):
        from keras.metrics import binary_crossentropy
        return binary_crossentropy

    def getOptimizer(self):
        from keras.optimizers import rmsprop
        return rmsprop()
