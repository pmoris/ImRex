""" CNN model for recognizing generated peptides. """
# import keras
from keras.models import Sequential
from keras.layers import (
    #    Dense,
    #    Dropout,
    #    Flatten,
    Conv2D,
    #    MaxPool2D,
    #    LeakyReLU,
    GlobalAveragePooling2D,
    Activation,
)
from keras.layers.normalization import BatchNormalization

from src.models.model import Model


NUM_CLASSES = 1


class ModelDebug(Model):
    def __init__(self, width, height, channels, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height
        self.channels = channels

    def _build_model(self):
        model = Sequential()

        input_shape = (self.width, self.height, self.channels)
        # WEIGHT_DECAY = 1e-6
        KERNEL_INIT = "he_normal"

        def create_conv(
            depth,
            kernel_size=(3, 3),
            activation="relu",
            padding="same",
            kernel_initializer=KERNEL_INIT,
            **kwargs
        ):
            return Conv2D(
                depth,
                kernel_size,
                activation=activation,
                padding=padding,
                kernel_initializer=kernel_initializer,
                **kwargs
            )

        model.add(create_conv(8, kernel_size=(3, 3), input_shape=input_shape))
        # model.add(Dropout(0.4))
        model.add(BatchNormalization())
        model.add(create_conv(8, kernel_size=(3, 3)))
        # model.add(Dropout(0.25))
        model.add(BatchNormalization())

        # model.add(create_conv(8, kernel_size=(3, 3)))
        # model.add(Dropout(0.4))
        # model.add(BatchNormalization())
        # model.add(create_conv(1, kernel_size=(1, 1)))
        # model.add(Dropout(0.25))
        # model.add(BatchNormalization())

        model.add(create_conv(NUM_CLASSES, (1, 1)))
        model.add(BatchNormalization())

        model.add(GlobalAveragePooling2D())
        model.add(Activation("sigmoid"))

        # model.add(Flatten())
        # model.add(Dense(32, activation='tanh'))
        # model.add(Dense(NUM_CLASSES, activation='sigmoid'))

        return model

    def get_loss(self):
        from keras.metrics import binary_crossentropy

        return binary_crossentropy

    def get_optimizer(self):
        from keras.optimizers import rmsprop

        return rmsprop()
