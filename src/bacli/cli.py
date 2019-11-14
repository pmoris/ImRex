import atexit
import inspect
import argparse


functions = dict()
description = ""


def setDescription(desc):
    global description
    description = desc


# decorator for runnable functions
def command(f):
    fName = f.__name__
    if fName in functions:
        print("Function {} already registered, overwriting it..")

    functions[fName] = f
    return f


def parse_args():
    parser = argparse.ArgumentParser(description=description)
    subparsers = parser.add_subparsers(help="Select one of the following subcommands:", dest='command', metavar="subcommand")
    subparsers.required = True

    for fName, f in functions.items():
        sub_parser = subparsers.add_parser(fName, help=f.__doc__, description=f.__doc__)
        for param in inspect.signature(f).parameters.values():
            tpe = param.annotation
            nargs = None
            if tpe is inspect.Parameter.empty:
                tpe = str
            elif tpe is list:
                nargs = "*"
                tpe = str
            if param.default is not inspect.Parameter.empty:
                prefix = "-" if len(param.name) == 1 else "--"
                sub_parser.add_argument(prefix + param.name,
                                        help="type: {}, default={}".format(tpe.__name__, param.default),
                                        type=tpe, default=param.default, nargs=nargs)
            else:
                sub_parser.add_argument(param.name, help="type: " + tpe.__name__,
                                        type=tpe, nargs=nargs)

    cmd_args = parser.parse_args()
    fName = cmd_args.command
    f = functions[fName]
    args = cmd_args._get_args()
    kwargs = {n: v for n, v in cmd_args._get_kwargs() if n != "command"}
    f(*args, **kwargs)


def main():
    try:
        if len(functions) > 0:
            parse_args()
    except argparse.ArgumentError:
        pass
    except SystemExit:
        pass


# if imported
if __name__ != "__main__":
    atexit.register(main)
