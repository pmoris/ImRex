""" CNN model for recognizing generated peptides. """
# import keras
from keras.layers import Conv2D, Dense, Dropout, Flatten, MaxPool2D  # , LeakyReLU
from keras.layers.normalization import BatchNormalization
from keras.models import Sequential

from src.models.model import Model


class ModelPPIPadded(Model):
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
            padding="valid",
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

        model.add(create_conv(16, kernel_size=(10, 10), input_shape=input_shape))
        # model.add(BatchNormalization())
        # model.add(create_conv(16, kernel_size=(5, 5)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        # model.add(create_conv(32, kernel_size=(5, 5), input_shape=input_shape))
        # model.add(BatchNormalization())
        model.add(create_conv(16, kernel_size=(8, 8)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        # model.add(create_conv(32, kernel_size=(3, 3), input_shape=input_shape))
        # model.add(BatchNormalization())
        model.add(create_conv(16, kernel_size=(5, 5)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        # model.add(create_conv(32, kernel_size=(3, 3), input_shape=input_shape))
        # model.add(BatchNormalization())
        model.add(create_conv(8, kernel_size=(3, 3)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        model.add(create_conv(8, kernel_size=(3, 3)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        model.add(Flatten())
        model.add(Dense(16, activation="tanh"))
        model.add(Dense(1, activation="sigmoid"))

        return model

    def get_loss(self):
        from keras.metrics import binary_crossentropy

        return binary_crossentropy

    def get_optimizer(self):
        from keras.optimizers import rmsprop

        return rmsprop()
