import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from QtWidgets.SampleListWidget import SampleListWidget
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


class CentralWidget(QWidget):
    def __init__(self, parent=None):
        super(CentralWidget, self).__init__(parent)

        self.setLayout(QGridLayout())

        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)
        self.axes.grid(True)
        self.axes.plot(range(1, 5), range(1, 5))

        self.canvas = FigureCanvas(self.fig)

        self.layout().addWidget(self.canvas)

        #self.plot_widget.axes.plot([0, 1, 2, 3, 4], [10, 1, 20, 3, 40])


class AppMainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(AppMainWindow, self).__init__(parent)
        self.setWindowTitle('Dock')
        self.resize(1400, 700)

        self.setLayout(QHBoxLayout())

        bar = self.menuBar()
        file_menu = bar.addMenu('File')
        file_menu.addAction('New')
        file_menu.addAction('Load')

        statusBar = QStatusBar()
        #p = statusBar.palette()
        #p.setColor(statusBar.backgroundRole(), Qt.red)
        #statusBar.setPalette(p)
        statusBar.setAttribute(Qt.WA_StyledBackground, True)
        statusBar.setStyleSheet('background-color: #939393;')
        self.setStatusBar(statusBar)

        self.file_list = SampleListWidget()
        self.file_list_dockable = QDockWidget('List Dockable', self)
        self.file_list_dockable.setWidget(self.file_list)

        self.file_list_dockable.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.file_list_dockable.setFloating(False)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.file_list_dockable)

        self.setCentralWidget(CentralWidget())