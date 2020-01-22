""" CNN model for recognizing generated peptides. """
# import keras
# import keras_metrics
# from keras.layers import Dense, Dropout, Flatten, Conv2D, MaxPool2D
# from keras.layers import LeakyReLU
# import keras.initializers
from keras import Sequential
from keras.layers import Dropout, Conv2D, MaxPool2D
from keras.layers.normalization import BatchNormalization
from keras.layers import Activation
from keras.layers import GlobalAveragePooling2D  # , GlobalMaxPooling2D
from keras.regularizers import l2

from src.models.model import Model

NUM_CLASSES = 1
LENGTH = 10


class ModelGAP(Model):
    def __init__(self, channels, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.channels = channels

    def _build_model(self):
        WEIGHT_DECAY = 1e-6
        # KERNEL_INIT = keras.initializers.he_normal

        KERNEL_INIT = "he_normal"

        def create_conv(
            depth,
            kernel_size=(3, 3),
            activation="relu",
            padding="same",
            regularizer=l2(WEIGHT_DECAY),
            kernel_initializer=KERNEL_INIT,
            **kwargs
        ):
            return Conv2D(
                depth,
                kernel_size,
                activation=activation,
                padding=padding,
                kernel_regularizer=regularizer,
                kernel_initializer=kernel_initializer,
                **kwargs
            )

        model = Sequential()

        model.add(create_conv(256, (3, 3), input_shape=(None, None, self.channels)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())
        model.add(create_conv(256, (3, 3)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())
        model.add(create_conv(128, (3, 3)))
        model.add(MaxPool2D(pool_size=(2, 2), padding="same"))
        # model.add(Dropout(0.25))
        model.add(BatchNormalization())

        # model.add(create_conv(64, (3, 3)))
        # model.add(BatchNormalization())
        # model.add(create_conv(32, (3, 3)))
        # model.add(BatchNormalization())
        # model.add(create_conv(16, (3, 3)))
        # # model.add(MaxPool2D(pool_size=(3, 3), strides=(2, 2), padding='same'))
        # model.add(MaxPool2D(pool_size=(2, 2), padding='same'))
        # model.add(Dropout(0.25))
        # model.add(BatchNormalization())

        model.add(create_conv(NUM_CLASSES, (1, 1)))
        model.add(BatchNormalization())

        model.add(GlobalAveragePooling2D())
        model.add(Activation("sigmoid"))

        return model

    def get_loss(self):
        from keras.metrics import binary_crossentropy

        return binary_crossentropy

    def get_optimizer(self):
        from keras.optimizers import rmsprop

        return rmsprop()
