""" CNN model for recognizing generated peptides. """
from typing import Optional

# from keras import Model
# from keras.layers import LeakyReLU
# from keras.layers import Activation
# from keras.regularizers import l2
import tensorflow
from tensorflow.keras.layers import (
    BatchNormalization,
    Conv1D,
    Dense,
    Dropout,
    GlobalMaxPooling1D,
    Input,
)
from tensorflow.keras.regularizers import l2

from src.models.model import Model

NUM_CLASSES = 1
LENGTH = 10


class ModelSeparatedInputsNetTcrSelu(Model):
    def __init__(
        self,
        optimizer,
        learning_rate: Optional[float] = None,
        regularization: Optional[float] = None,
        dropout_conv: Optional[float] = None,
        dropout_dense: Optional[float] = None,
        activation_function: str = "selu",
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.optimizer = optimizer
        self.learning_rate = learning_rate
        self.regularization = l2(regularization) if regularization else None
        self.dropout_conv = dropout_conv
        self.dropout_dense = dropout_dense
        self.activation_function = activation_function.lower()

    def _build_model(self):
        if self.activation_function == "selu":
            KERNEL_INIT = "lecun_normal"
        else:
            KERNEL_INIT = "he_normal"

        input1 = Input(shape=(None, 20))
        input2 = Input(shape=(None, 20))

        def feature_extraction(input):  # noqa: A002
            convolutions = list()
            convolutions.append(
                Conv1D(
                    100,
                    (1,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT,
                    activation=self.activation_function,
                    kernel_regularizer=self.regularization,
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (3,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT,
                    activation=self.activation_function,
                    kernel_regularizer=self.regularization,
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (5,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT,
                    activation=self.activation_function,
                    kernel_regularizer=self.regularization,
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (7,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT,
                    activation=self.activation_function,
                    kernel_regularizer=self.regularization,
                )(input)
            )
            convolutions.append(
                Conv1D(
                    100,
                    (9,),
                    padding="same",
                    kernel_initializer=KERNEL_INIT,
                    activation=self.activation_function,
                    kernel_regularizer=self.regularization,
                )(input)
            )

            return convolutions

        conv_input1 = feature_extraction(input1)
        conv_input2 = feature_extraction(input2)

        # second convolutional layer - concat pairwise length-wise per kernel size
        concat_conv_list = [
            tensorflow.keras.layers.concatenate([conv_1, conv_2], axis=1)
            for conv_1, conv_2 in zip(conv_input1, conv_input2)
        ]

        batchnorm_1 = [BatchNormalization()(concat) for concat in concat_conv_list]

        conv_2 = [
            tensorflow.keras.layers.Conv1D(
                filters=100,
                kernel_size=1,
                kernel_initializer=KERNEL_INIT,
                padding="same",
                activation=self.activation_function,
                kernel_regularizer=self.regularization,
            )(bn)
            for bn in batchnorm_1
        ]

        max_pool = [GlobalMaxPooling1D()(conv) for conv in conv_2]

        if self.dropout_conv:
            drop_out_conv = [Dropout(self.dropout_conv)(mp) for mp in max_pool]
            batch_norm_2 = [BatchNormalization()(d) for d in drop_out_conv]
        else:
            batch_norm_2 = [BatchNormalization()(mp) for mp in max_pool]

        concat_2 = tensorflow.keras.layers.concatenate(batch_norm_2, axis=1)

        dense = Dense(
            10,
            activation=self.activation_function,
            kernel_regularizer=self.regularization,
        )(concat_2)
        if self.dropout_dense:
            drop_out_dense = Dropout(self.dropout_dense)(dense)
            predictions = Dense(1, activation="sigmoid")(drop_out_dense)
        else:
            predictions = Dense(1, activation="sigmoid")(dense)

        model = tensorflow.keras.Model(inputs=[input1, input2], outputs=predictions)

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
