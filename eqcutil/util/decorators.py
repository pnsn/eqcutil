import os, csv
from pathlib import Path

def catch_args(func, *args, **kwargs):
    callargs = inspect.getcallargs(func, *args, **kwargs)
    callargs.pop('self')
    _kwargs = callargs.pop("kwargs", {})
    info = {'args': callargs,
            'kwargs': _kwargs}
    result = func(*args, **kwargs)
    self = args[0]
    self._internal_catch(info)
    return result

def save_kwargs(mode='w'):
    """Decorator that saves the key-word arguments and key-word arguments
    of a decorated function to the specified file

    Modified from a prompt output to ChatGPT

    :param mode: write mode, defaults to 'w'
    :type mode: str, optional
    """
    def decorator(func):
        def wrapper(**kwargs):
            arg_names = func.__code__.co_varnames[:func.__code__.co_argcount]

            row = [kwargs.get(name, None) for name in arg_names[len(kwargs):]]
            filename = Path().cwd() / f'{func.__name__}_kwargs.csv'

            with open(filename, mode) as file:
                writer = csv.writer(file)
                writer.writerow(arg_names)
                writer.writerow(row)
            return func(**kwargs)
        return wrapper
    return decorator