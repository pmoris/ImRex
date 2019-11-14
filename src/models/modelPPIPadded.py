""" CNN model for recognizing generated peptides. """
from .model import Model


class ModelPPIPadded(Model):
    def __init__(self, width, height, channels, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height
        self.channels = channels

    def _buildModel(self):
        import keras
        from keras.models import Sequential
        from keras.layers import Dense, Dropout, Flatten, Conv2D, MaxPool2D, LeakyReLU
        from keras.layers.normalization import BatchNormalization

        model = Sequential()

        inputShape = (self.width, self.height, self.channels)
        WEIGHT_DECAY = 1e-6
        KERNEL_INIT = "he_normal"

        def createConv(
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

        model.add(createConv(16, kernel_size=(10, 10), input_shape=inputShape))
        # model.add(BatchNormalization())
        # model.add(createConv(16, kernel_size=(5, 5)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        # model.add(createConv(32, kernel_size=(5, 5), input_shape=inputShape))
        # model.add(BatchNormalization())
        model.add(createConv(16, kernel_size=(8, 8)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        # model.add(createConv(32, kernel_size=(3, 3), input_shape=inputShape))
        # model.add(BatchNormalization())
        model.add(createConv(16, kernel_size=(5, 5)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        # model.add(createConv(32, kernel_size=(3, 3), input_shape=inputShape))
        # model.add(BatchNormalization())
        model.add(createConv(8, kernel_size=(3, 3)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        model.add(createConv(8, kernel_size=(3, 3)))
        model.add(MaxPool2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))
        model.add(BatchNormalization())

        model.add(Flatten())
        model.add(Dense(16, activation="tanh"))
        model.add(Dense(1, activation="sigmoid"))

        return model

    def getLoss(self):
        from keras.metrics import binary_crossentropy

        return binary_crossentropy

    def getOptimizer(self):
        from keras.optimizers import rmsprop

        return rmsprop()
