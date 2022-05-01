#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Launch the Concatenator GUI"""

import multiprocessing
from itaxotools.concatenator_gui import run

if __name__ == '__main__':
    multiprocessing.freeze_support()
    run()
