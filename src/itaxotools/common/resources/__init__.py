# This must be a package for importlib.resources and pyinstaller
# to detect and import the contained data files.

import pathlib
import sys
import os

__all__ = ['package_path', 'get_common', 'get_local', 'get']


def package_path_pyinstaller(package):
    if isinstance(package, str):
        parts = package.split('.')
    elif isinstance(package, type(sys)):
        parts = package.__name__.split('.')
    else:
        return None
    path = pathlib.Path(sys._MEIPASS)
    for part in parts:
        path /= part
    return path


def package_path_file(package):
    if isinstance(package, str):
        file = sys.modules[package].__file__
    elif isinstance(package, type(sys)):
        file = package.__file__
    else:
        return None
    path = pathlib.Path(os.path.dirname(file))
    return path


def package_path_importlib(package):
    return importlib.resources.files(package)


try:
    import importlib.resources2
    package_path = package_path_importlib
except ModuleNotFoundError:
    if hasattr(sys, '_MEIPASS'):
        package_path = package_path_pyinstaller
    else:
        package_path = package_path_file

_resource_path = package_path(__package__)


def get_common(path):
    return str(_resource_path / path)


def get_local(package, path):
    return str(package_path(package) / path)


def get(*args):
    if len(args) == 1:
        return get_common(*args)
    return get_local(*args)
