import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from QtWidgets.MainWindow import AppMainWindow


if __name__ == '__main__':
    app = QApplication(sys.argv)
    MainWindow = AppMainWindow()
    MainWindow.show()
    sys.exit(app.exec_())