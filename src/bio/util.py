from pathlib import Path

import numpy as np


# decorator annotation
def decorator(func):
    """ Decorator for decorator functions.  # noqa: D401.

    Allows for syntax:

    @decorator
    def foo(func, kwarg1=None):
        pass    # noqa: RST3-1

    which works both when a parameter list is given when using the decorator foo or not.

        @foo
        def bar():
            pass    # noqa: RST3-1

        @foo(5)
        def bar():
            pass    # noqa: RST3-1

    """

    def wrapper(*args, **kwargs):
        # If no parameters given (only function)
        if len(args) == 1 and len(kwargs) == 0 and callable(args[0]):
            return func(*args, **kwargs)

        # If parameters given
        def decorator(f):
            return func(f, *args, **kwargs)

        return decorator

    return wrapper


def after(postprocess):
    """ Debug decorator. """

    def wrapper(f):
        def wrapped(*args, **kwargs):
            res = f(*args, **kwargs)
            postprocess(res)
            return res

        return wrapped

    return wrapper


def scale_matrix(
    matrix: np.ndarray,
    oldlower: float,
    oldupper: float,
    lower: float = 0.0,
    upper: float = 255.0,
) -> np.ndarray:
    """Scale matrix inplace to new bounds, based on a provided initial minimum and maximum value.

    E.g. a matrix with combined pairwise amino acid features could be rescaled,
    while using the min and max combined value for all possible combinations among the existing 20 amino acids,
    rather than just the min and max that are present in the given matrix.

    Note: the matrix will be modified in-place, as well as returning it.
    See: https://stackoverflow.com/questions/10149416/numpy-modify-array-in-place
    This is caused by the in-place operators. Statements like matrix = matrix + x, would not cause this behaviour.

    Note: In-place operations do not change the dtype of the container array.
    Since the desired scaled values are floats,
    the matrix needs to have the correct dtype before the in-place operations are performed.

    Note: sklearn.preprocessing.MinMaxScaler or optimized code like `*= 255.0/image.max()`
    (source: https://stackoverflow.com/questions/1735025/how-to-normalize-a-numpy-array-to-within-a-certain-range)
    is not used, because in these implementations the old minimum and maximum values
    cannot be specified, but are automatically derived from the supplied data.

    Parameters
    ----------
    matrix : np.ndarray
        A matrix with pairwise amino acid properties. Provided by bio.peptide_feature.matrix
        or bio.operator.matrix.
    oldlower : float
        Minimum theoretical value that the amino acid property can take,
        provided by bio.operator.Operator().min_op().
    oldupper : float
        Maximum theoretical value that the amino acid property can take,
        provided by bio.operator.Operator().max_op().
    lower : float
        Minimum value to use for rescaling, default = 0.
    upper : float
        Maximum value to use for rescaling.

    Returns
    -------
    np.ndarray
        The matrix with all the values scaled between the new lower and upper bounds,
        given the supplied theoretical lower and upper bounds (e.g. based on the
        min and max combination of any amino acid for the given operator).
    """
    matrix -= oldlower
    matrix *= upper - lower
    matrix /= oldupper - oldlower
    matrix += lower
    return matrix


def subdirs(directory):
    """Return all subdirectories of a given directory path."""
    return [x for x in Path(directory).iterdir() if x.is_dir()]
