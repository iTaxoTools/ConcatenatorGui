# PyInstaller entry points for setuptools
# https://pyinstaller.readthedocs.io/en/stable/hooks.html

# When a module is detected by PyInstaller, it will search
# for corresponding hooks and tests in this directory.

import os

def get_hook_dirs():
    return [os.path.dirname(__file__)]

def get_PyInstaller_tests():
    return [os.path.dirname(__file__)]
