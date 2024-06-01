#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt6 import QtWidgets, uic
from PyQt6.QtCore import QTimer, Qt
from PyQt6.QtGui import QPainter, QBrush, QPen
import pyqtgraph as pg


import sys
import numpy as np
import os
import logging

from .funcs import *

pg.setConfigOption('background', 'k')
pg.setConfigOption('foreground', 'w')

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

file_dir = os.path.dirname(__file__)

main_win_ui = os.path.join(file_dir, 'PorousModel.ui')
param_widget_ui = os.path.join(file_dir, 'paramWidget.ui')

class ParamWidget(QtWidgets.QWidget):
    def __init__(self, parent, callback):
        super().__init__(parent)
        uic.loadUi(param_widget_ui, self)
        
        self.callback = callback
        self.sigmaBox.setMaximum(10000)
        self.sigmaBox.setValue(15)
        self.thicknessBox.setMaximum(10000)
        self.thicknessBox.setValue(10)
        self.phiBox.setMaximum(5.0)
        self.phiBox.setValue(0.99)
        self.tortuosityBox.setMaximum(5.0)
        self.tortuosityBox.setValue(1.01)
        self.vLengthBox.setMaximum(100000)
        self.vLengthBox.setValue(87.4)
        self.tLengthBox.setMaximum(100000)
        self.tLengthBox.setValue(196)
        
    def get_alpha(self):
        sigma = self.sigmaBox.value() * 1e3
        d = self.thicknessBox.value() * 1e-3
        tort = self.tortuosityBox.value()
        phi = self.phiBox.value()
        v_length = self.vLengthBox.value() * 1e-6
        t_length = self.tLengthBox.value() * 1e-6
        
        if self.DBMRadio.isChecked():
            alpha = DBM(freq, d, sigma)
        elif self.JCARadio.isChecked():
            alpha = JCA(freq, d,  phi, sigma, tort, v_length, t_length)
        else:
            alpha = DBM(freq, d, sigma)
            
        return alpha
    
    def updatePlot(self):
        self.callback()
        

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__() # Call the inherited classes __init__ method
        uic.loadUi(main_win_ui, self) # Load the .ui file
        self.color = 0
        
        
        self.setup()
        
        self.show() # Show the GUI
        
    def addModelCurve(self):
        widget = ParamWidget(self, self.updatePlot)
        self.paramWidgets.append(widget)
        i = len(self.paramWidgets)
        self.tabWidget.addTab(widget, "Model {}".format(i))
        self.absorptionCurves.append(self.ax.plot(freq, np.zeros_like(freq), pen=self.color, name='Model {}'.format(i)))
        self.color+=1 
        self.updatePlot()
        
    def setup(self):
        self.ax = self.graphWidget.addPlot(title="Absorption")        
        self.paramWidgets = []
        self.absorptionCurves = []
        self.addModelCurve()
        self.ax.enableAutoRange('xy', False)
        self.ax.setYRange(0, 1, padding=0.1)
        self.ax.showGrid(x=True, y=True)
        
    def updatePlot(self):
        for i, widget in enumerate(self.paramWidgets):
            alpha = widget.get_alpha()
                
            self.absorptionCurves[i].setData(freq, alpha)
        
def main():
    logging.basicConfig()
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    app.exec()
        
if __name__ == '__main__':
    main()