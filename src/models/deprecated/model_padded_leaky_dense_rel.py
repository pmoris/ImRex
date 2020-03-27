""" CNN model for recognizing generated peptides. """
# import keras
from tensorflow.keras.layers import (
    BatchNormalization,
    Conv2D,
    Dense,
    Dropout,
    Flatten,
    LeakyReLU,
    MaxPool2D,
)
from tensorflow.keras.models import Sequential

from src.models.model import Model


NUM_CLASSES = 1


class ModelPaddedLeakyDenseRel(Model):
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
            # activation="lrelu",
            padding="same",
            kernel_initializer=KERNEL_INIT,
            **kwargs
        ):
            return Conv2D(
                depth,
                kernel_size,
                # activation=activation,
                padding=padding,
                kernel_initializer=kernel_initializer,
                **kwargs
            )

        model.add(create_conv(128, kernel_size=(3, 3), input_shape=input_shape))
        model.add(LeakyReLU(alpha=0.3))
        # model.add(Dropout(0.4))
        model.add(BatchNormalization())
        model.add(create_conv(64, kernel_size=(3, 3)))
        model.add(LeakyReLU(alpha=0.3))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        # model.add(Dropout(0.4))
        model.add(BatchNormalization())

        model.add(create_conv(128, kernel_size=(3, 3)))
        model.add(LeakyReLU(alpha=0.3))
        # model.add(Dropout(0.4))
        model.add(BatchNormalization())
        model.add(create_conv(64, kernel_size=(3, 3)))
        model.add(LeakyReLU(alpha=0.3))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        # model.add(Dropout(0.4))
        model.add(BatchNormalization())

        model.add(Flatten())
        model.add(Dense(32))
        model.add(LeakyReLU(alpha=0.3))
        model.add(Dense(NUM_CLASSES, activation="sigmoid"))

        return model

    def get_loss(self):
        from tensorflow.keras.losses import BinaryCrossentropy

        return BinaryCrossentropy()

    def get_optimizer(self):
        from tensorflow.keras.optimizers import RMSprop

        return RMSprop()
