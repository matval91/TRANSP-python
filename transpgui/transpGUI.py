#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:08:00 2019

@author: vallar
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication
import sys
import transpGUI_design as design
import pytransp.classes.transp_exp as te
import numpy as np


class TranspGUI(QtWidgets.QMainWindow, design.Ui_MainWindow):
    def __init__(self, parent=None):
        super(TranspGUI, self).__init__(parent)
        self.setupUi(self)
        
        self.ReadTranspButton.clicked.connect(self.ReadTranspButtonPressed)
        self.listWidget1d.itemClicked.connect(self.listWidget1dClicked)
        
        
    def ReadTranspButtonPressed(self):
        options = QtWidgets.QFileDialog.Options()
        self.fname, _ = QtWidgets.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","CDF files(*.CDF)", options=options)        
        _translate = QtCore.QCoreApplication.translate
        self.FileNameLabel.setText(_translate("MainWindow", self.fname))
        self.transp = te.transp_exp(self.fname)
        self._updatelists()

    def _updatelists(self):
        variables=self.transp.file.variables.keys()
        self.dict_1d = {}
        self.dict_2d = {}
        for v in np.sort(variables):
            tmpv = self.transp.file.variables[v]
            n_dims = np.size(np.shape(tmpv))
            #item = QtWidgets.QListWidgetItem()
            #item.setData(v)
            #item.setText('%s: %s' % (v,self.transp.file.variables[v].long_name))
            if n_dims == 1:
                self.listWidget1d.addItem('%s: %s' % (v,self.transp.file.variables[v].long_name))
                self.dict_1d[v] = self.transp.file.variables[v].long_name
                #self.listWidget1d.addItem(item)
            elif n_dims ==2:
                self.listWidget2d.addItem('%s %s' % (v,self.transp.file.variables[v].long_name))
                self.dict_2d[v] = self.transp.file.variables[v].long_name
                #self.listWidget2d.addItem(item)
            else:
                continue
    
    def listWidget1dClicked(self):
        item = self.listWidget1d.currentItem()
        
    
def main():
    app = QApplication(sys.argv)
    form = TranspGUI()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()