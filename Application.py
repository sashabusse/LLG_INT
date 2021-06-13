import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from QtWidgets.MainWindow import AppMainWindow


class Application(QApplication):
    def __init__(self, ):
        super(Application, self).__init__(sys.argv)
        