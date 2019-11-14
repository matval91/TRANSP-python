# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'transpGUI_design.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1189, 684)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.Variables = QtWidgets.QToolBox(self.centralwidget)
        self.Variables.setGeometry(QtCore.QRect(10, 90, 481, 331))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.Variables.setFont(font)
        self.Variables.setObjectName("Variables")
        self.pagelist1d = QtWidgets.QWidget()
        self.pagelist1d.setGeometry(QtCore.QRect(0, 0, 481, 257))
        self.pagelist1d.setObjectName("pagelist1d")
        self.listWidget1d = QtWidgets.QListWidget(self.pagelist1d)
        self.listWidget1d.setGeometry(QtCore.QRect(10, -10, 451, 271))
        self.listWidget1d.setObjectName("listWidget1d")
        self.Variables.addItem(self.pagelist1d, "")
        self.pagelist2d = QtWidgets.QWidget()
        self.pagelist2d.setGeometry(QtCore.QRect(0, 0, 471, 297))
        self.pagelist2d.setObjectName("pagelist2d")
        self.listWidget2d = QtWidgets.QListWidget(self.pagelist2d)
        self.listWidget2d.setGeometry(QtCore.QRect(10, -10, 441, 271))
        self.listWidget2d.setObjectName("listWidget2d")
        self.Variables.addItem(self.pagelist2d, "")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setGeometry(QtCore.QRect(0, 0, 681, 80))
        font = QtGui.QFont()
        font.setPointSize(13)
        self.groupBox.setFont(font)
        self.groupBox.setObjectName("groupBox")
        self.ReadTranspButton = QtWidgets.QPushButton(self.groupBox)
        self.ReadTranspButton.setGeometry(QtCore.QRect(10, 40, 171, 31))
        self.ReadTranspButton.setObjectName("ReadTranspButton")
        self.FileNameLabel = QtWidgets.QLabel(self.groupBox)
        self.FileNameLabel.setGeometry(QtCore.QRect(230, 40, 441, 31))
        self.FileNameLabel.setObjectName("FileNameLabel")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1189, 20))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.Variables.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.Variables.setItemText(self.Variables.indexOf(self.pagelist1d), _translate("MainWindow", "1D"))
        self.Variables.setItemText(self.Variables.indexOf(self.pagelist2d), _translate("MainWindow", "2D"))
        self.groupBox.setTitle(_translate("MainWindow", "Read file"))
        self.ReadTranspButton.setText(_translate("MainWindow", "Read CDF file"))
        self.FileNameLabel.setText(_translate("MainWindow", "First, Read Transp File "))

