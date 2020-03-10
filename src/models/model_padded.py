""" CNN model for recognizing generated peptides. """
from typing import Optional

from keras.layers import Conv2D, Dense, Dropout, Flatten, MaxPool2D  # , LeakyReLU
from keras.layers.normalization import BatchNormalization
from keras.models import Sequential

from src.models.model import Model


NUM_CLASSES = 1


class ModelPadded(Model):
    def __init__(
        self,
        width: int,
        height: int,
        channels: int,
        optimizer: str,
        learning_rate: Optional[bool] = None,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height
        self.channels = channels
        self.optimizer = optimizer
        self.learning_rate = learning_rate

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
            **kwargs,
        ):
            return Conv2D(
                depth,
                kernel_size,
                activation=activation,
                padding=padding,
                kernel_initializer=kernel_initializer,
                **kwargs,
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
        model.add(Dense(32, activation="tanh"))  # relu?
        model.add(Dense(NUM_CLASSES, activation="sigmoid"))

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
