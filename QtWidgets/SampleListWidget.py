import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *


class SampleNode:
    def __init__(self, sample, name):
        self.sample = sample
        self.name = name


class SampleListWidget(QWidget):
    def __init__(self, parent=None):
        super(SampleListWidget, self).__init__(parent)
        self.setLayout(QVBoxLayout())

        self.list_widget = QListWidget()

        self.layout().addWidget(self.list_widget)

    def add_node(self, node):
        item = QListWidgetItem()
        item.setData(Qt.UserRole, node)
        item.setData(Qt.DisplayRole, node.name)

        self.list_widget.addItem(item)

    def remove_selected(self):
        items = self.list_widget.selectedItems()
        for item in items:
            self.list_widget.takeItem(self.list_widget.row(item))