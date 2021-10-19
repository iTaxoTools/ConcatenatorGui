# https://pyinstaller.readthedocs.io/en/stable/hooks.html

from PyInstaller.utils.hooks import collect_data_files
datas = collect_data_files('itaxotools.concatenator_gui')
