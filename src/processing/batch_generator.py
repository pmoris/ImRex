from keras.utils import Sequence
import numpy as np

from .stream import BatchStream


class BatchGenerator(Sequence, BatchStream):
    def __init__(self, batchStream, batchSize, multipleInput=False):
        BatchStream.__init__(self, batchStream)
        Sequence.__init__(self)
        self.batchStream = batchStream
        self.batchSize = batchSize
        self.multipleInput = multipleInput

    def getBatch(self, batchSize, *args, **kwargs):
        # print("Batch requested")
        batch = self.batchStream.getBatch(batchSize, *args, **kwargs)
        # transform from    (X1, y1), (X2, y2), ...
        #           to      (X1, X2, ...), (y1, y2, ...)
        # print(batch)
        X, y = list(zip(*batch))
        # print("Batch delivered")
        if self.multipleInput:
            X = list(zip(*X))
            ret = [np.array(x) for x in X], np.array(y)
            # print(ret)
            return ret
        else:
            return np.array(X), np.array(y)

    def __getitem__(self, index):
        """Gets batch at position `index`.

        # Arguments
            index: position of the batch in the Sequence.

        # Returns
            A batch
        """
        # print(f"\nStart batch {index}")
        batch = self.getBatch(self.batchSize)
        # print(f"\nDone batch {index}")
        return batch

    def __len__(self):
        return len(self.batchStream) // self.batchSize

    def on_epoch_end(self):
        """Method called at the end of every epoch.
        """
        self.sendEvent("epoch_end")
