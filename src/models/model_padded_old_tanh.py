""" CNN model for recognizing generated peptides. """
from tensorflow.keras.layers import (
    BatchNormalization,
    Conv2D,
    Dense,
    Dropout,
    Flatten,
    MaxPool2D,
)  # , LeakyReLU
from tensorflow.keras.models import Sequential

from src.models.model import Model


NUM_CLASSES = 1


class ModelPadded(Model):
    def __init__(self, width, height, channels, optimizer, include_lr, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height
        self.channels = channels
        self.optimizer = optimizer
        self.include_lr = include_lr

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

        model.add(create_conv(128, kernel_size=(3, 3), input_shape=input_shape))
        # model.add(Dropout(0.4))
        model.add(BatchNormalization())
        model.add(create_conv(64, kernel_size=(3, 3)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        # model.add(Dropout(0.4))
        model.add(BatchNormalization())

        model.add(create_conv(128, kernel_size=(3, 3)))
        # model.add(Dropout(0.4))
        model.add(BatchNormalization())
        model.add(create_conv(64, kernel_size=(3, 3)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        # model.add(Dropout(0.4))
        model.add(BatchNormalization())

        model.add(Flatten())
        model.add(Dense(32, activation="tanh"))
        model.add(Dense(NUM_CLASSES, activation="sigmoid"))

        return model

    def get_loss(self):
        from tensorflow.keras.losses import BinaryCrossentropy

        return BinaryCrossentropy()

    def get_optimizer(self):
        if self.optimizer == "rmsprop":

            from tensorflow.keras.optimizers import RMSprop

            return RMSprop()

        elif self.optimizer == "adam":

            from tensorflow.keras.optimizers import Adam

            if self.include_lr:
                return Adam(lr=self.include_lr)
            else:
                return Adam()

        elif self.optimizer == "SGD":

            from tensorflow.keras.optimizers import SGD

            if self.include_lr:
                return SGD(lr=self.include_lr)
            else:
                return SGD()
