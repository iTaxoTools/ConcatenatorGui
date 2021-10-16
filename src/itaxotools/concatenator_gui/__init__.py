
"""GUI entry point"""


def run():
    """
    Show the Concatenator window and enter the main event loop.
    Imports are made locally to optimize multiprocessing.
    """

    import sys
    from PySide6 import QtWidgets
    from PySide6 import QtCore
    from .main import Main

    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('Fusion')
    files = (file for file in sys.argv[1:])
    main = Main(files=files)
    main.setWindowFlags(QtCore.Qt.Window)
    main.setModal(True)
    main.show()
    sys.exit(app.exec_())
