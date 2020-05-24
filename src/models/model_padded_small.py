""" CNN model for recognizing generated peptides. """
from typing import Optional

from tensorflow.keras.layers import (
    BatchNormalization,
    Conv2D,
    Dense,
    Dropout,
    Flatten,
    GlobalAveragePooling2D,
    MaxPool2D,
    SpatialDropout2D,
)
from tensorflow.keras.models import Sequential
from tensorflow.keras.regularizers import l2

from src.models.model import Model


NUM_CLASSES = 1


class ModelPaddedSmall(Model):
    def __init__(
        self,
        width: int,
        height: int,
        channels: int,
        optimizer: str,
        depth1: int,
        depth2: int,
        gap: bool = False,
        learning_rate: Optional[float] = None,
        regularization: Optional[float] = None,
        activation_function_conv: str = "relu",
        activation_function_dense: str = "tanh",
        dropout_conv: Optional[float] = None,
        dropout_dense: Optional[float] = None,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height
        self.channels = channels
        self.optimizer = optimizer.lower()
        self.learning_rate = learning_rate
        self.depth1 = depth1
        self.depth2 = depth2
        self.gap = gap
        self.regularization = l2(regularization) if regularization else None
        self.activation_function_conv = activation_function_conv.lower()
        self.activation_function_dense = activation_function_dense.lower()
        self.dropout_conv = dropout_conv
        self.dropout_dense = dropout_dense

    def _build_model(self):
        model = Sequential()

        input_shape = (self.width, self.height, self.channels)
        # WEIGHT_DECAY = 1e-6

        if self.activation_function_conv == "selu":
            KERNEL_INIT = "lecun_normal"
        else:
            KERNEL_INIT = "he_normal"

        def create_conv(
            depth,
            kernel_size=(3, 3),
            activation=self.activation_function_conv,
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
                kernel_regularizer=self.regularization,
                **kwargs,
            )

        model.add(create_conv(self.depth1, kernel_size=(3, 3), input_shape=input_shape))
        model.add(BatchNormalization())
        model.add(create_conv(self.depth1, kernel_size=(3, 3)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        if self.dropout_conv:
            model.add(SpatialDropout2D(self.dropout_conv))
        model.add(BatchNormalization())

        model.add(create_conv(self.depth2, kernel_size=(3, 3)))
        model.add(BatchNormalization())
        model.add(create_conv(self.depth2, kernel_size=(3, 3)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        if self.dropout_conv:
            model.add(SpatialDropout2D(self.dropout_conv))
        model.add(BatchNormalization())

        if self.gap:
            model.add(GlobalAveragePooling2D())
        else:
            model.add(Flatten())

        model.add(
            Dense(
                32,
                activation=self.activation_function_dense,
                kernel_regularizer=self.regularization,
            )
        )
        if self.dropout_dense:
            model.add(Dropout(self.dropout_dense))

        model.add(Dense(NUM_CLASSES, activation="sigmoid"))

        return model

    def get_loss(self):
        from tensorflow.keras.losses import BinaryCrossentropy

        return BinaryCrossentropy()

    def get_optimizer(self):
        if self.optimizer == "rmsprop":

            from tensorflow.keras.optimizers import RMSprop

            if self.learning_rate:
                return RMSprop(learning_rate=self.learning_rate)
            else:
                return RMSprop()

        elif self.optimizer == "adam":

            from tensorflow.keras.optimizers import Adam

            if self.learning_rate:
                return Adam(learning_rate=self.learning_rate)
            else:
                return Adam()

        elif self.optimizer == "SGD":

            from tensorflow.keras.optimizers import SGD

            if self.learning_rate:
                return SGD(llearning_rate=self.learning_rate)
            else:
                return SGD()
