from typing import List

from src.processing.data_stream import DataStream


def tee(stream: DataStream, amount: int = 2) -> List[DataStream]:
    """ Replicates a given stream amount times.

    Note that this will consume the provided DataStream, or
    in the case of for example Filter, will consume the underlying
    DataStream object. This will cause get() or next() calls to
    raise a StopIteration error, but will not prevent loops or lists
    from being made.

    Parameters
    ----------
    stream : DataStream
        The DataStream to be replicated.
    amount : int, optional
        The number of times it should be replicated, by default 2

    Returns
    -------
    list
        A list of replicated DataStreams.
    """
    data = list(stream)
    return [DataStream(data) for _ in range(amount)]
