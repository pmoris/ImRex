import os


# decorator annotation
def decorator(func):
    """
    Decorator for decorator functions. Allows for syntax:
    ```
    @decorator
    def foo(func, kwarg1=None):
        pass
    ```
    which works both when a parameter list is given when using the decorator foo or not.
    ```
    @foo
    def bar():
        pass

    @foo(5)
    def bar():
        pass

    ```

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


def scale_matrix(matrix, oldlower, oldupper, upper=255.0):
    """ Scales matrix inplace to new bounds """
    # matrix = matrix.astype(float)
    lower = 0.0
    matrix -= oldlower  # - lower
    matrix *= (upper - lower) / (oldupper - oldlower)
    # matrix += lower
    return matrix


def subdirs(directory):
    return filter(
        lambda x: os.path.isdir(x),
        [os.path.join(directory, subdir) for subdir in os.listdir(directory)],
    )
