""" CNN model for recognizing generated peptides. """
from .model import Model

import keras_metrics



class ModelPPILit(Model):
    def __init__(self, width, height, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.height = height

    def _buildModel(self):
        import keras
        from keras import Model, Input
        from keras.layers import Dense, Dropout, Flatten, Conv1D, MaxPool1D, Embedding, LSTM, Activation
        from keras.layers.normalization import BatchNormalization
        from keras.regularizers import l2
        import keras.initializers

        WEIGHT_DECAY = 1e-6
        KERNEL_INIT = keras.initializers.he_normal

        input1 = Input(shape=(self.width,))
        input2 = Input(shape=(self.height,))

        def featureExtraction(input):
            embedding = Embedding(21, 128)(input)
            x = Conv1D(10, 10, padding='valid',kernel_initializer=KERNEL_INIT())(embedding)
            x = Activation("relu")(x)
            x = MaxPool1D(pool_size=2, padding='valid')(x)
            x = Conv1D(10, 8, padding='valid', kernel_initializer=KERNEL_INIT())(x)
            x = Activation("relu")(x)
            x = MaxPool1D(pool_size=2, padding='valid')(x)
            x = Conv1D(10, 5, padding='valid', kernel_initializer=KERNEL_INIT())(x)
            x = Activation("relu")(x)
            x = MaxPool1D(pool_size=2, padding='valid')(x)
            out = LSTM(80)(x)
            # out = Flatten()(x)
            return out

        part1 = featureExtraction(input1)
        part2 = featureExtraction(input2)

        merged_vector = keras.layers.concatenate([part1, part2], axis=-1)
        predictions = Dense(1, activation='sigmoid')(merged_vector)

        model = Model(inputs=[input1, input2], outputs=predictions)
        return model

    def getLoss(self):
        from keras.metrics import binary_crossentropy
        return binary_crossentropy

    def getOptimizer(self):
        from keras.optimizers import adam
        return adam()
