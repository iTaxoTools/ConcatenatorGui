
"""GUI entry point"""


def main():
    """
    Show the Concatenator main window.
    Imports are made locally to optimize multiprocessing.
    """

    import sys
    import pathlib
    from PySide6 import QtWidgets
    from PySide6 import QtCore
    from .main import Main

    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('Fusion')
    main = Main()
    main.setWindowFlags(QtCore.Qt.Window)
    main.setModal(True)
    main.show()
    sys.exit(app.exec_())
