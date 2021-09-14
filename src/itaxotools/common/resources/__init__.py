# This must be a package for importlib.resources and pyinstaller
# to detect and import the contained data files.

import pathlib
import sys

try:
    import importlib.resources
    _resource_path = importlib.resources.files(__package__)
except ModuleNotFoundError as e:
    if hasattr(sys, '_MEIPASS'):
        _resource_path = (pathlib.Path(sys._MEIPASS) /
            'itaxotools' / 'common' / 'resources')
    else:
        import os
        _resource_path = pathlib.Path(os.path.dirname(__file__))

def get(path):
    return str(_resource_path / path)

def get_icon(path):
    return str(_resource_path / 'icons/svg' / path)
