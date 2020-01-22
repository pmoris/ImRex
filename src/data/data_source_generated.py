import numpy as np

# import src.bacli as bacli
# from src.neural.trainer import createCheckpointer, createLRR, createTensorboardCallback


# BATCH_SIZE = 256
# EPOCHS = 20
#
# # settings for generator
# TRAIN_STEPS = 50
# VALIDATION_STEPS = 20


class GeneratedDataSource(object):
    def generate_sample(self):
        raise NotImplementedError

    def generate_batch(self, batch_size):
        X = list()
        Y = list()
        for i in range(batch_size):
            x, y = self.generate_sample()
            X.append(x)
            Y.append(y)

        return np.array(X), np.array(Y).reshape((batch_size, 1))

    def generate_data(self, batch_size=BATCH_SIZE):
        while True:
            yield self.generate_batch(batch_size)


# def createCommands(dataSource):
#     @bacli.command
#     def train(modelFile: str, epochs: int = 40, batch_size: int = BATCH_SIZE,
#               trainSteps=TRAIN_STEPS, validationSteps=VALIDATION_STEPS,
#               includeLRR:bool = True):
#         """ Train the CNN. """
#         from keras.models import load_model
#
#         model = load_model(modelFile)
#
#         print("Training model:")
#         model.summary()
#
#         callbacks = list()
#         callbacks.append(createCheckpointer(modelFile))
#         callbacks.append(createTensorboardCallback(os.path.splitext(modelFile)[0]))
#         if includeLRR:
#             callbacks.append(createLRR())
#
#         print("Fitting CNN")
#         h = model.fit_generator(generator=dataSource.generateData(batch_size),
#                                 steps_per_epoch=trainSteps,
#                                 epochs=epochs,
#                                 validation_data=dataSource.generateData(batch_size),
#                                 validation_steps=validationSteps,
#                                 verbose=1,
#                                 callbacks=callbacks)
