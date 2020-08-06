import numpy as np
from tensorflow.keras.utils import Sequence

from src.processing.stream import BatchStream


class BatchGenerator(Sequence, BatchStream):
    def __init__(self, batch_stream, batch_size, multiple_input=False):
        BatchStream.__init__(self, batch_stream)
        Sequence.__init__(self)
        self.batch_stream = batch_stream
        self.batch_size = batch_size
        self.multiple_input = multiple_input

    def get_batch(self, batch_size, *args, **kwargs):
        # print("Batch requested")
        batch = self.batch_stream.get_batch(batch_size, *args, **kwargs)
        # transform from    (X1, y1), (X2, y2), ...
        #           to      (X1, X2, ...), (y1, y2, ...)
        # print(batch)
        X, y = list(zip(*batch))
        # print("Batch delivered")
        if self.multiple_input:
            X = list(zip(*X))
            ret = [np.array(x) for x in X], np.array(y)
            # print(ret)
            return ret
        else:
            return np.array(X), np.array(y)

    def __getitem__(self, index):
        """Get batch at position `index`.

        # Arguments
            index: position of the batch in the Sequence.

        # Returns
            A batch
        """
        # print(f"\nStart batch {index}")
        batch = self.get_batch(self.batch_size)
        # print(f"\nDone batch {index}")
        return batch

    def __len__(self):
        return len(self.batch_stream) // self.batch_size

    def on_epoch_end(self):
        """Method called at the end of every epoch."""  # noqa: D401
        self.send_event("epoch_end")
