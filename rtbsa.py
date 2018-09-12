#!/usr/local/lcls/package/python/current/bin/python
#Written by Zimmer, modified by Lisa

import sys
import os
from matplotlib.backends.backend_qt4 import FigureCanvasQT as FigureCanvas
from epics import caget, ca, PV
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.figure import Figure
from xml.etree import ElementTree
from re import sub
from shutil import copy
from time import sleep
import numpy as np
from numpy import polyfit,poly1d,polyval,corrcoef,std,delete,mean
from math import floor
from scipy.optimize import curve_fit
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from subprocess import Popen, PIPE, check_output
import pyqtgraph as pg
import pyqtgraph.exporters

from rtbsa_ui import Ui_RTBSA as BSA_UI

class RTBSA(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.ui = BSA_UI()
        self.ui.setupUi(self)
        self.setWindowTitle('Real Time BSA')
        self.loadStyleSheet()
        self.setUpGraph()        
        QObject.connect(self.ui.enter1, SIGNAL("textChanged(const QString&)"),
                        self.search1)
        QObject.connect(self.ui.enter2, SIGNAL("textChanged(const QString&)"),
                        self.search2)
        self.ui.listWidget.itemClicked.connect(self.setenter1)
        self.ui.listWidget_2.itemClicked.connect(self.setenter2)
        self.bsapvs = ['GDET:FEE1:241:ENRC', 'GDET:FEE1:242:ENRC',
                       'GDET:FEE1:361:ENRC', 'GDET:FEE1:362:ENRC']
        
        #Generate list of BSA PVS
        bsapvs = check_output(['eget','-s','LCLS:BSA_PVS']).splitlines()[1:-1]
        
        self.bsapvs.extend([pv.split()[-2] for pv in bsapvs])
        for pv in self.bsapvs:
            self.ui.listWidget.addItem(pv)
            self.ui.listWidget_2.addItem(pv)
        
        commonlist = ['GDET:FEE1:241:ENRC', 'GDET:FEE1:242:ENRC',
                      'GDET:FEE1:361:ENRC', 'GDET:FEE1:362:ENRC',
                      'KLYS:LI20:K6:VOLT', 'ACCL:IN20:300:L0A_P',
                      'ACCL:IN20:300:L0A_A', 'KLYS:LI20:K7:VOLT',
                      'ACCL:IN20:400:L0B_P', 'ACCL:IN20:400:L0B_A',
                      'KLYS:LI20:K8:VOLT', 'BPMS:IN20:731:X',
                      'ACCL:LI21:1:L1S_P', 'ACCL:LI21:1:L1S_A',
                      'KLYS:LI21:K1:VOLT', 'ACCL:LI21:180:L1X_P',
                      'ACCL:LI21:180:L1X_A', 'KLYS:LI21:K2:VOLT',
                      'BPMS:LI21:233:X', 'BLEN:LI21:265:AIMAX',
                      'BPMS:LI24:801:X', 'BLEN:LI24:886:BIMAX',
                      'BPMS:LTU1:250:X', 'BPMS:LTU1:450:X', 'BPMS:UND1:1090:X',
                      'BPMS:UND1:1090:Y', 'BPMS:UND1:2090:X',
                      'BPMS:UND1:2090:Y', 'BLD:SYS0:500:UND_POS_X',
                      'BLD:SYS0:500:UND_ANG_X', 'BLD:SYS0:500:UND_POS_Y',
                      'BLD:SYS0:500:UND_ANG_Y']
                      
        self.ui.common1.addItems(commonlist)
        self.ui.common2.addItems(commonlist)
        self.ui.common1.setCurrentIndex(24)
        self.ui.common1.activated.connect(self.commonactivated)
        self.ui.common2.activated.connect(self.commonactivated)
        self.ui.AvsB.clicked.connect(self.AvsBClick)
        self.ui.draw_button.clicked.connect(self.on_draw)        
        self.ui.stop_button.clicked.connect(self.stop)
        self.ui.log_button.clicked.connect(self.logbook)
        self.ui.mcclog_button.clicked.connect(self.MCCLog)
        self.ui.avg_cb.clicked.connect(self.avg_click)
        self.ui.std_cb.clicked.connect(self.std_click)
        self.ui.corr_cb.clicked.connect(self.corr_click)
        self.ui.parab_cb.clicked.connect(self.parab_click)
        self.ui.line_cb.clicked.connect(self.line_click)
        self.ui.fitedit.returnPressed.connect(self.fitorderactivated)       
        self.ui.common1_rb.clicked.connect(self.common_1_click)
        self.ui.common2_rb.clicked.connect(self.common_2_click)
        self.ui.enter1_rb.clicked.connect(self.enter_1_click)
        self.ui.enter2_rb.clicked.connect(self.enter_2_click)
        self.ui.AvsT_cb.clicked.connect(self.AvsTClick)
        self.ui.AFFT.clicked.connect(self.AFFTClick)
        self.ui.enter1.returnPressed.connect(self.commonactivated)
        self.ui.enter2.returnPressed.connect(self.commonactivated)
        self.ui.points.returnPressed.connect(self.points_entered)
        
        #Initial number of points
        self.numpoints=2800
        #20ms polling time
        self.updatetime=20
        #Set initial polynomial fit to 2
        self.fitorder=2
        self.dpi = 100
        
        self.ui.fitedit.setDisabled(True)
        self.ui.enter1.setDisabled(True)
        self.ui.enter2.setDisabled(True)
        self.ui.label.setDisabled(True)
        self.ui.listWidget.setDisabled(True)
        self.ui.listWidget_2.setDisabled(True)
        self.statusBar().showMessage('Hi there!  I missed you!')
        self.abort=True
        self.ui.parab_cb.setChecked(False)
        #Used to update plot
        self.timer=QTimer(self)
        self.rate = PV('EVNT:SYS0:1:LCLSBEAMRATE') 
        #Kludge to prevent double fit plotting
        #self.timer2=QTimer(self) 
        self.menuBar().setStyleSheet('QWidget{background-color:grey;color:purple}') 
        self.create_menu()
        self.create_status_bar()

    def setUpGraph(self):
        self.plot = pg.PlotWidget(alpha=0.75)
        #self.p1 = self.plot.plotItem
        layout = QGridLayout()
        self.ui.widget.setLayout(layout)
        layout.addWidget(self.plot,0,0)    
        self.plot.showGrid(1,1)
        #self.pg_item = self.plot.getPlotItem()
        #self.curve = pg.PlotCurveItem(pen='o')
        #self.curve=self.plot.plot()
        #self.plot.addItem(self.curve)

    def loadStyleSheet(self):
        try:
            self.cssfile = "/home/physics/zimmerc/python/style.css"
            with open(self.cssfile,"r") as f:
                self.setStyleSheet(f.read())
        except:
            print "Error loading style sheet"
            pass

    def create_status_bar(self):
        self.status_text = QLabel()
        palette = QPalette()
        palette.setColor(palette.Foreground,Qt.magenta)
        self.statusBar().addWidget(self.status_text, 1)
        self.statusBar().setPalette(palette)        

    def search(self, enter, widget):
        widget.clear()
        query = str(enter.text())
        for pv in self.bsapvs:
            if query.lower() in pv.lower():
                widget.addItem(pv)
    
    def search1(self):
        search(self.ui.enter1, self.ui.listWidget)
        
    def search2(self):
        search(self.ui.enter2, self.ui.listWidget_2)

    def setEnter(self, widget, enter, search, enter_rb):
        selection = widget.currentItem()
        enter.textChanged.disconnect()
        enter.setText(selection.text())
        QObject.connect(enter, SIGNAL("textChanged(const QString&)"), search)
        if not self.abort and enter_rb.isChecked():
            self.stop()
            self.on_draw()
    
    def setenter1(self):
        setEnter(self.ui.listWidget, self.ui.enter1, self.search1,
                 self.ui.enter1_rb)

    def setenter2(self):
        setEnter(self.ui.listWidget_2, self.ui.enter2, self.search2,
                 self.ui.enter2_rb)

    def points_entered(self):
        try:
            self.numpoints = int(self.ui.points.text())
        except ValueError:
            self.statusBar().showMessage('Enter an integer, 1 to 2800', 6000)
            self.numpoints = 120
            self.ui.points.setText('120')
            return
        
        if self.numpoints > 2800:
            self.statusBar().showMessage('Max # points is 2800', 6000)
            self.numpoints = 2800
            self.ui.points.setText('2800')
            return
            
        if self.numpoints < 1:
            self.statusBar().showMessage('Min # points is 1', 6000)
            self.numpoints = 1
            self.ui.points.setText('1')
            return
            
        self.reinitialize_plot()

    def pv1select(self):
        self.pv1 = self.ui.pv1box.currentText()
        if self.pv1 == 'Select PV':
            return
        else:
            self.stop()
            
    def appendToPvList(self, common_rb, common, enter_rb, enter, device):
        success = True
        
        if common_rb.isChecked():
            self.pvlist.append(str(common.currentText()+'HSTBR'))
            
        elif enter_rb.isChecked():
            if str(enter.text()).strip():
                self.pvlist.append(str(enter.text())+'HSTBR')
            else:
                self.statusBar().showMessage('Device '+ device 
                                             + ' field blank. Aborting.', 10000)
                self.ui.draw_button.setEnabled(True)
                success = False
            
        return success
        
    def setValSynced(self):
    
        numBadShots = round((self.time2 - self.time1)*self.rate.value)
        
        val1synced = self.val1pre[max(0, numBadShots)
                                  :min(2800, 2800 + numBadShots)]
        val2synced = self.val2pre[max(0, -numBadShots)
                                  :min(2800, 2800 - numBadShots)]
            
        return [abs(numBadShots), val1synced, val2synced]
        
    def updateRate(self):
        rate = self.rate.value
        while rate not in [120.0, 60.0, 30.0, 10.0]:
            QApplication.processEvents()
            rate = self.rate.value
        return rate
 
    def cleanPlot(self):
        self.plot.clear()
        
        self.avg_text = pg.TextItem('', color=(200, 200, 250), anchor=(0, 1))
        self.std_text = pg.TextItem('', color=(200, 200, 250), anchor=(0, 1))
        self.slope_text = pg.TextItem('', color=(200, 200, 250), anchor=(0, 1))
        self.corr_text = pg.TextItem('', color=(200, 200, 250), anchor=(0, 1))

        plotLabels = [self.avg_text, self.std_text, self.slope_text,
                      self.corr_text]
             
        for plotLabel in plotLabels:
            self.plot.addItem(plotLabel)
            
    def initialzeData(self):
        self.statusBar().showMessage('Initializing...')
            
        if self.ui.common1_rb.isChecked():
            self.device = str(self.ui.common1.currentText() + 'HSTBR')
            
        elif self.ui.enter1_rb.isChecked():
            if str(self.ui.enter1.text()).strip():
                self.device = str(self.ui.enter1.text() + 'HSTBR')
        else:
            #self.statusBar().showMessage('Device A field blank. Aborting.',
                                          #10000)
            #self.ui.draw_button.setEnabled(True)
            return None
        
        # If i use Popen instead of pyepics caget, the program doesn't
        # start lagging if you change PVs (?!?!?). Some stupid bug in 
        # new pyepics.
        #newdata=caget(self.device)
        getdata = Popen("caget " + self.device, stdout=PIPE, shell=True)
        newdata = str(getdata.communicate()).split()[2:-1]
        newdata[-1] = newdata[-1][:-4]
        return [float(i) for i in newdata]
        
    def plotFit(self, newdata):
        self.curve = pg.PlotCurveItem(newdata[2800-self.numpoints:2800], pen=1)
        self.plot.addItem(self.curve)
        self.plot.setTitle(self.device)
        self.xdata=range(self.numpoints)
        
        #Fit line
        if self.ui.line_cb.isChecked():
            (m, b) = polyfit(self.xdata, newdata[2800 - self.numpoints:2800], 1)
            fitdata = polyval([m, b], self.xdata)
            m = "{:.3e}".format(m)
            self.slope_text.setText('Slope: ' + str(m))    
            self.fit = pg.PlotCurveItem(self.xdata,fitdata, 'g-', linewidth = 1)
            self.plot.addItem(self.fit)
        
        #Fit polynomial
        elif self.ui.parab_cb.isChecked():
            self.ui.fitedit.setDisabled(False)
            co = polyfit(self.xdata, newdata[2800-self.numpoints:2800],
                         self.fitorder)
            pol = poly1d(co)
            fit = pol(self.xdata)
            self.parab = pg.PlotCurveItem(self.xdata, fit, pen=3)
            self.plot.addItem(self.parab)
            if self.fitorder == 2:
                self.slope_text.setText('Peak: '+str(-co[1]/(2*co[0])))
            elif self.fitorder == 3:
                self.slope_text.setText(str("{:.2e}".format(co[0])) + 'x^3' 
                      + str("+{:.2e}".format(co[1])) + 'x^2' 
                      + str("+{:.2e}".format(co[2])) + 'x' 
                      + str("+{:.2e}".format(co[3])))  
                
    def genTimePlotA(self):
        newdata = self.initialzeData()
            
        # Using pyepics caget causes progressive slowdown of GUI!?!?
        #newdata=caget(self.device)
        
        if not newdata:
            self.statusBar().showMessage('Invalid PV? Unable to get data.'
                       + ' Aborting.', 10000)
            self.ui.draw_button.setEnabled(True)
            return
        try:
            self.plotFit(newdata)
        except UnboundLocalError:
            self.statusBar().showMessage('No Data, Aborting Plotting Algorithm',
                                         10000)
            return
            
        self.statusBar().showMessage('Running')
        self.timer=QTimer(self)
        self.timer.singleShot(self.updatetime, self.update_plot_HSTBR)
        
    def adjustVals(self):
        self.updateRate()
        numBadShots, val1synced, val2synced = self.setValSynced()
            
        blength = 2800 - numBadShots
        if (self.numpoints <= blength):
            self.val1 = val1synced[blength-self.numpoints:blength]     
            self.val2 = val2synced[blength-self.numpoints:blength]
        else:
            self.val1=val1synced
            self.val2=val2synced
        
    def updateValsFromInput(self):
        self.pvlist = []
            
        if not self.appendToPvList(self.ui.common1_rb, self.ui.common1,
                                   self.ui.enter1_rb, self.ui.enter1, "A"):
            return
                          
        if not self.appendToPvList(self.ui.common2_rb, self.ui.common2,
                                   self.ui.enter2_rb, self.ui.enter2, "B"):
            return
        
        self.val1, self.val2 = [], []
        
        self.val1pv = PV(self.pvlist[0], form='time')
        self.val2pv = PV(self.pvlist[1], form='time')
        
        self.statusBar().showMessage('Initializing/Syncing (be patient, '
                                     + 'may take 5 seconds)...')
        
        self.time1, self.time2 = 1, 2
        
        # Are there other rates besides 0 and 1?
        
            
        self.booyah = self.val1pv.add_callback(self.ItDoneChanged) 
        self.booyah2 = self.val2pv.add_callback(self.ItDoneChanged2)
        
        while (self.time1 == 1 or self.time2 == 2) and not self.abort:
            QApplication.processEvents()
            
        self.adjustVals()
            
    
    # TODO This is near identical to plotFit - I'll ask Chris if the differences
    # are necessary
    def genABPlot(self):
        self.curve = pg.ScatterPlotItem(self.val1, self.val2, pen=1, symbol='x',
                                        size=5)
        self.plot.addItem(self.curve)
        self.plot.setTitle(self.pvlist[1] + ' vs. ' + self.pvlist[0])
        
        if self.ui.line_cb.isChecked():
            self.fit = pg.PlotCurveItem(pen=3)
            try:
                (m, b) = polyfit(self.val1, self.val2, 1)
                fitdata = polyval([m, b], self.val1)
                m = "{:.3e}".format(m)
                self.slope_text.setText('Slope: ' + str(m))    
                self.fit.setData(self.val1, fitdata)
            except:
                print "Error getting line fit"
                pass
            self.plot.addItem(self.fit)
            
        elif self.ui.parab_cb.isChecked():
            self.parab = pg.PlotCurveItem(pen=3, size=2)
            try:
                co = polyfit(self.val1, self.val2, self.fitorder)
                pol = poly1d(co)
                sorted1 = sorted(self.val1)
                fit = pol(sorted1)
                self.parab.setData(sorted1, fit)
                
                if self.fitorder == 2:
                    self.slope_text.setText('Peak: ' + str(-co[1]/(2*co[0])))
                    
                elif self.fitorder == 3:
                    self.slope_text.setText(str("{:.2e}".format(co[0])) + 'x^3'
                          + str("+{:.2e}".format(co[1])) + 'x^2'
                          + str("+{:.2e}".format(co[2])) + 'x' 
                          + str("+{:.2e}".format(co[3])))
                
            except np.linalg.linalg.LinAlgError:
                print "Linear algebra error getting curve fit"
                pass
            except:
                self.slope_text.setText('Fit failed')
                pass
            self.plot.addItem(self.parab)
    
    def genFFTPlot(self):
        newdata = self.initialzeData()
        newdata = newdata[2800-self.numpoints:2800]
        newdata.extend(np.zeros(self.numpoints*2).tolist())
        newdata = newdata - np.mean(newdata);
        #ps = np.log(np.abs(np.fft.fft(newdata[2800-self.numpoints:2800]))
                    #/self.numpoints)
        ps = np.abs(np.fft.fft(newdata))/len(newdata)
        self.FS = self.rate.value
        self.freqs = np.fft.fftfreq(len(newdata), 1.0/self.FS)
        self.keep = self.freqs >= 0
        ps = ps[self.keep]
        self.freqs = self.freqs[self.keep]
        self.idx = np.argsort(self.freqs)
        self.curve = pg.PlotCurveItem(x = self.freqs[self.idx],
                                      y= ps[self.idx], pen=1)
        self.plot.addItem(self.curve)
        self.plot.setTitle(self.device)
        #self.plot.setTitle(self.pvlist[1]+' vs. '+self.pvlist[0])
    
    def genPlotAndSetTimer(self, genPlot):
        if self.abort:
            return False
            
        try:
            genPlot()
        except UnboundLocalError:
            self.statusBar().showMessage('No Data, Aborting Plotting Algorithm',
                                         10000)
            return False
    
        self.timer=QTimer(self)
        self.timer.singleShot(self.updatetime, self.update_BSA_Plot)
        self.statusBar().showMessage('Running')
        return True
    
    # Where the magic happens(well, where it starts to happen). This initializes
    # the BSA plotting and then starts a timer to update the plot.
    def on_draw(self):
        plotTypeIsValid = (self.ui.AvsT_cb.isChecked() or self.ui.AvsB.isChecked()
                           or self.ui.AFFT.isChecked()) 
        
        if not plotTypeIsValid:
            self.statusBar().showMessage('Pick a Plot Type (PV vs. time or B vs A)',
                                         10000)
            return
            
        self.ui.draw_button.setDisabled(True)
        self.abort = False
        
        self.cleanPlot()
            
        ####Plot history buffer for one PV####
        if self.ui.AvsT_cb.isChecked():
            self.genTimePlotA()
        
        ####Plot for 2 PVs####
        elif self.ui.AvsB.isChecked(): 
            self.updateValsFromInput()
            
            if not self.genPlotAndSetTimer(self.genABPlot):
                return
        
        ####Plot power spectrum####
        else: 
            if not self.genPlotAndSetTimer(self.genFFTPlot):
                return
    
    def update_BSA_Plot(self):
        QApplication.processEvents()
        
        if self.abort:
            return
            
        self.plot.showGrid(self.ui.grid_cb.isChecked(),
                           self.ui.grid_cb.isChecked())
        
        self.adjustVals()
        
        #Filter out NaNs and ridiculous BLEN values
        elements_to_use = np.ones(len(self.val2), dtype=bool)
        elements_to_use = np.logical_and(elements_to_use, 
                 np.logical_not(np.isnan(self.val2)))
        elements_to_use = np.logical_and(elements_to_use,
                 np.logical_not(np.isnan(self.val1)))
                 
        if "BLEN:LI24:886" in self.pvlist[1]:
            elements_to_use = np.logical_and(elements_to_use, 
                     np.logical_not(self.val2 > 12000))
        if "BLEN:LI24:886" in self.pvlist[0]:
            elements_to_use = np.logical_and(elements_to_use, 
                        np.logical_not(self.val1 > 12000))
                        
        self.val1 = self.val1[elements_to_use]
        self.val2 = self.val2[elements_to_use]
        
        self.curve.setData(self.val1, self.val2)
        
        if self.ui.autoscale_cb.isChecked():
            mx = np.max(self.val2)
            mn = np.min(self.val2)
            
            if mn != mx:
                self.plot.setYRange(mn, mx)
                
            mx = np.max(self.val1)
            mn = np.min(self.val1)
            
            if mn != mx:
                self.plot.setXRange(mn, mx)
                
        if self.ui.avg_cb.isChecked():
            a = mean(self.val2)
            a = "{:.4}".format(a)
            self.avg_text.setPos(min(self.val1), min(self.val2))
            self.avg_text.setText('AVG: ' + str(a))
            
        if self.ui.std_cb.isChecked():
            s = std(self.val2)
            s = "{:.4}".format(s)
            #TODO Is this a typo?
            self.std_text.setPos((min(self.val1) + (min(self.val1) 
                   + max(self.val1))/2)/2,
                  min(self.val2))
            self.std_text.setText('STD: ' + str(s))
            
        if self.ui.corr_cb.isChecked():
            correlation = corrcoef(self.val1, self.val2)
            self.corr_text.setPos(min(self.val1), max(self.val2))
            self.corr_text.setText("Corr. Coefficient: {:.3}"
                    .format(correlation.item(1)))
                    
        if self.ui.line_cb.isChecked():
            (m, b) = polyfit(self.val1, self.val2, 1)   
            fitdata = polyval([m, b], self.val1)
            m = "{:.4}".format(m)
            self.slope_text.setPos((min(self.val1) + max(self.val1))/2,
                    min(self.val2))
            self.slope_text.setText('Slope: ' + str(m))        
            self.fit.setData(self.val1, fitdata)
            
        elif self.ui.parab_cb.isChecked():
            try:
                co = polyfit(self.val1, self.val2, self.fitorder)
                pol = poly1d(co)
                fit = pol(sorted(self.val1))
                self.parab.setData(sorted(self.val1), fit)
                
                if self.fitorder == 2:
                    peak = -co[1]/(2*co[0])
                    peak = "{:.4}".format(peak)
                    self.slope_text.setText('Peak: ' + str(peak))
                    
                elif self.fitorder == 3:
                    self.slope_text.setText(str("{:.2e}".format(co[0])) + 'x^3'
                          + str("+{:.2e}".format(co[1])) + 'x^2'
                          + str("+{:.2e}".format(co[2])) + 'x'
                          + str("+{:.2e}".format(co[3])))
                          
                self.slope_text.setPos((min(self.val1) + max(self.val1))/2,
                        min(self.val2))
            except:
                print "Error getting fit"
                pass
        self.timer.singleShot(self.updatetime, self.update_BSA_Plot)

    def update_plot_HSTBR(self):
        
        self.plot.showGrid(self.ui.grid_cb.isChecked(),
                           self.ui.grid_cb.isChecked())
        
        QApplication.processEvents()
        
        if self.abort:
            return
            
        #newdata=caget(self.device)
        getdata = Popen("caget " + self.device, stdout=PIPE, shell=True)
        newdata = str(getdata.communicate()).split()[2:-1]
        newdata[-1] = newdata[-1][:-4]
        newdata = [float(i) for i in newdata]
        chopped = newdata[2800 - self.numpoints:2800]
        self.curve.setData(chopped)
        
        if self.ui.autoscale_cb.isChecked():
            mx = max(chopped)
            mn = min(chopped)
            if mx - mn > .00001:
                self.plot.setYRange(mn, mx)
                self.plot.setXRange(0, len(chopped))
                
        if self.ui.avg_cb.isChecked():
            a = mean(chopped)
            a = "{:.3e}".format(a)
            self.avg_text.setPos(0, min(chopped))
            self.avg_text.setText('AVG: ' + str(a))
            
        if self.ui.std_cb.isChecked():
            s = std(chopped)
            s = "{:.3e}".format(s)   
            self.std_text.setPos(self.numpoints/4, min(chopped))
            self.std_text.setText('STD: ' + str(s))
            
        if self.ui.corr_cb.isChecked():
            self.corr_text.setText('')
            
        if self.ui.line_cb.isChecked():
            (m, b) = polyfit(self.xdata, chopped, 1)
            fitdata=polyval([m, b], self.xdata)
            m = "{:.3e}".format(m)
            self.slope_text.setPos(self.numpoints/2, min(chopped))
            self.slope_text.setText('Slope: ' + str(m))
            self.fit.setData(self.xdata, fitdata)
            
        elif self.ui.parab_cb.isChecked():
            try:
                co = polyfit(self.xdata, chopped, self.fitorder)
                pol = poly1d(co)
                fit = pol(self.xdata)
                self.parab.setData(self.xdata, fit)
                
                if self.fitorder == 2:
                    peak = -co[1]/(2*co[0])
                    peak = "{:.3e}".format(peak)
                    self.slope_text.setText('Peak: ' + str(peak))
                    
                elif self.fitorder == 3:
                    self.slope_text.setText(str("{:.2e}".format(co[0])) + 'x^3'
                          + str("+{:.2e}".format(co[1])) + 'x^2'
                          + str("+{:.2e}".format(co[2])) + 'x'
                          + str("+{:.2e}".format(co[3])))
                self.slope_text.setPos(self.numpoints/2,min(chopped))
            
            except:
                print "Error getting fit"
                pass
                
        self.timer.singleShot(40, self.update_plot_HSTBR)

    def update_plot_FFT(self):
        self.plot.showGrid(self.ui.grid_cb.isChecked(), 
               self.ui.grid_cb.isChecked())
        QApplication.processEvents()
        
        if self.abort:
            return
            
        getdata = Popen("caget " + self.device,stdout=PIPE, shell=True)
        newdata = str(getdata.communicate()).split()[2:-1]
        newdata[-1] = newdata[-1][:-4]
        newdata = [float(i) for i in newdata]
        newdata = newdata[2800-self.numpoints:2800]
        newdata = np.array(newdata)
        nans, x = np.isnan(newdata), lambda z: z.nonzero()[0]
        #interpolate nans
        newdata[nans] = np.interp(x(nans), x(~nans), newdata[~nans]) 
        newdata = newdata - np.mean(newdata);
        newdata = newdata.tolist()
        newdata.extend(np.zeros(self.numpoints*2).tolist())
        ps = np.abs(np.fft.fft(newdata))/len(newdata)
        self.FS = self.rate.value
        self.freqs = np.fft.fftfreq(len(newdata), 1.0/self.FS)
        self.keep = (self.freqs >= 0)
        ps = ps[self.keep]
        self.freqs = self.freqs[self.keep]
        self.idx = np.argsort(self.freqs)
        self.curve.setData(x = self.freqs[self.idx], y = ps[self.idx])
        
        if self.ui.autoscale_cb.isChecked():
            mx = max(ps)
            mn = min(ps)
            if mx - mn >.00001:
                self.plot.setYRange(mn, mx)
                self.plot.setXRange(min(self.freqs), max(self.freqs))
                
        self.timer.singleShot(40, self.update_plot_FFT)
    
    #Callback function for PV1
    def ItDoneChanged(self, pvname=None, value=None, timestamp=None, **kw):
        self.time1 = timestamp
        self.val1pre = value    

    #Callback function for PV2
    def ItDoneChanged2(self, pvname=None, value=None, timestamp=None, **kw):
        self.time2 = timestamp
        self.val2pre = value

    def AvsTClick(self):
        palette = QPalette()
        if not self.ui.AvsT_cb.isChecked():
            pass
        else:
            self.ui.AvsB.setChecked(False)
            self.ui.AFFT.setChecked(False)
            self.AvsBClick()
        
    def AvsBClick(self):
        if not self.ui.AvsB.isChecked():
            self.ui.groupBox_2.setDisabled(True)
            self.ui.enter2_rb.setChecked(False)
            self.ui.enter2_rb.setDisabled(True)
            self.ui.enter2.setDisabled(True)
            self.ui.common2.setDisabled(True)
            self.ui.common2_rb.setChecked(False)
            self.ui.common2_rb.setDisabled(True)
        else:
            self.ui.AvsT_cb.setChecked(False)
            self.ui.AFFT.setChecked(False)
            self.AvsTClick()
            self.ui.groupBox_2.setDisabled(False)
            self.ui.listWidget_2.setDisabled(True)
            self.ui.enter2_rb.setDisabled(False)
            self.ui.enter2.setDisabled(True)
            self.ui.common2_rb.setDisabled(False)
            self.ui.common2_rb.setChecked(True)
            self.ui.common2.setDisabled(False)
        self.stop()

    def AFFTClick(self):
        if not self.ui.AFFT.isChecked():
            pass
        else:
            self.ui.AvsB.setChecked(False)
            self.ui.AvsT_cb.setChecked(False)
            self.AvsBClick()


    def avg_click(self):
        if not self.ui.avg_cb.isChecked():
           self.avg_text.setText('')

    def std_click(self):
        if not self.ui.std_cb.isChecked():
           self.std_text.setText('')

    def corr_click(self):
        if not self.ui.corr_cb.isChecked():
            self.corr_text.setText('')

    def enter_1_click(self):
        if self.ui.enter1_rb.isChecked():
            self.ui.enter1.setDisabled(False)
            self.ui.listWidget.setDisabled(False)
            self.ui.common1_rb.setChecked(False)
            self.ui.common1.setDisabled(True)
        else:
            self.ui.enter1.setDisabled(True)

    def enter_2_click(self):
        if self.ui.enter2_rb.isChecked():
            self.ui.enter2.setDisabled(False)
            self.ui.listWidget_2.setDisabled(False)
            self.ui.common2_rb.setChecked(False)
            self.ui.common2.setDisabled(True)
        else:
            self.ui.enter2.setDisabled(True)

    def common_1_click(self):
        if self.ui.common1_rb.isChecked():
            self.ui.common1.setEnabled(True)
            self.ui.enter1_rb.setChecked(False)
            self.ui.enter1.setDisabled(True)
            self.ui.listWidget.setDisabled(True)
        else:
            self.ui.common1.setEnabled(False)
        self.commonactivated()

    def commonactivated(self):
        if not self.abort:
            self.stop()
            self.timer.singleShot(150,self.on_draw)
            
            # Mechanism to fix a race condition- start a timer which allows the 
            # main loop to regain control and see that the abort variable has 
            # turned positive, and then the timer restarts the drawing process 
            # after a short time, giving the active loop time to stop.  This 
            # fixes a problem I was seeing (Gibbs noticed first) with a double 
            # plot showing up (rare occurence, only happened if the user tried 
            # to plot a blank PV name from the enter line)
            
            #QObject.connect(self.timer2, SIGNAL("timeout()"), self.on_draw)
            #self.timer2.start(self.updatetime)

    def common_2_click(self):
        if self.ui.common2_rb.isChecked():
            self.ui.common2.setEnabled(True)
            self.ui.enter2_rb.setChecked(False)
            self.ui.enter2.setDisabled(True)
            self.ui.listWidget_2.setDisabled(True)
        else:
            self.ui.common2.setEnabled(False)
        self.commonactivated()

    def line_click(self): 
        self.ui.parab_cb.setChecked(False)
        self.ui.fitedit.setDisabled(True)
        self.ui.label.setDisabled(True)
        self.reinitialize_plot()

    def fitorderactivated(self):
        try:
           self.fitorder = int(self.ui.fitedit.text())     
        except ValueError:
            self.statusBar().showMessage('Enter an integer, 1-30',6000)
            return
            
        if self.fitorder not in range(11) or self.fitorder == 0:
            self.statusBar().showMessage('Really?  That is going to be useful'
                    + ' to you?  The (already ridiculous)'
                    + ' range is 1-10.  Hope you win a '
                    + 'nobel prize jackass.', 6000)
            self.ui.fitedit.setText('2')
            self.fitorder = 2
            
        if self.fitorder != 2:
            try:
                self.slope_text.setText('')
            except AttributeError:
                pass
                        
    def parab_click(self):
        self.ui.line_cb.setChecked(False)
        
        if not self.ui.parab_cb.isChecked():
            self.ui.fitedit.setDisabled(True)
            self.ui.label.setDisabled(True)
        else:
            self.ui.fitedit.setEnabled(True)
            self.ui.label.setEnabled(True)
        self.reinitialize_plot()

    # This is a mess, but it works (used if user changes number points, 
    # fit type etc.) 
    def reinitialize_plot(self):
        self.cleanPlot()
        
        try:
         #Setup for single PV plotting 
            if self.ui.AvsT_cb.isChecked():
                self.genTimePlotA()
                    
            elif self.ui.AvsB.isChecked():
                self.genABPlot()
            else:
                getdata = Popen("caget " + self.device, stdout=PIPE, shell=True)
                newdata = str(getdata.communicate()).split()[2:-1]
                newdata[-1] = newdata[-1][:-4]
                newdata = [float(i) for i in newdata]
                newdata = newdata[2800 - self.numpoints:2800]
                newdata.extend(np.zeros(self.numpoints*2).tolist())
                newdata = newdata - np.mean(newdata);
                ps = np.abs(np.fft.fft(newdata))/len(newdata)
                rate = {4:1, 5:10, 6:30, 7: 60, 8:120}
                i = caget('IOC:BSY0:MP01:BYKIK_RATE')
                self.FS = rate[i];
                self.freqs = np.fft.fftfreq(self.numpoints, 1.0/self.FS)
                self.keep = (self.freqs >= 0)
                self.freqs = self.freqs[self.keep]
                ps = ps[self.keep]
                self.idx = np.argsort(self.freqs)
                self.curve = pg.PlotCurveItem(x = self.freqs[self.idx],
                         y = ps[self.idx], pen=1)
                self.plot.addItem(self.curve)
                self.plot.setTitle(self.device)

        except:
            pass

    # Shamelessly stolen from Shawn (thanks buddy). A lot of this is probably 
    # unnecessary for my purposes but I'm too lazy to clean it up
    def logbook(self):
        curr_time = datetime.now()
        timeString = curr_time.strftime("%Y-%m-%dT%H:%M:%S")
        log_entry = ElementTree.Element(None)
        severity  = ElementTree.SubElement(log_entry, 'severity')
        location  = ElementTree.SubElement(log_entry, 'location')
        keywords  = ElementTree.SubElement(log_entry, 'keywords')
        time      = ElementTree.SubElement(log_entry, 'time')
        isodate   = ElementTree.SubElement(log_entry, 'isodate')
        log_user  = ElementTree.SubElement(log_entry, 'author')
        category  = ElementTree.SubElement(log_entry, 'category')
        title     = ElementTree.SubElement(log_entry, 'title')
        metainfo  = ElementTree.SubElement(log_entry, 'metainfo')
        imageFile = ElementTree.SubElement(log_entry, 'link')
        imageFile.text = timeString + '-00.ps'
        thumbnail = ElementTree.SubElement(log_entry, 'file')
        thumbnail.text = timeString + "-00.png"
        text      = ElementTree.SubElement(log_entry, 'text')
        log_entry.attrib['type'] = "LOGENTRY"
        category.text = "USERLOG"
        location.text = "not set"
        severity.text = "NONE"
        keywords.text = "none"
        time.text = curr_time.strftime("%H:%M:%S")
        isodate.text =  curr_time.strftime("%Y-%m-%d")
        metainfo.text = timeString + "-00.xml"
        fileName = "/tmp/" + metainfo.text
        fileName = fileName.rstrip(".xml")
        log_user.text = 'Python Real-Time BSA'
        title.text = 'BSA Data'
        text.text = str(self.numpoints) + ' points'
        
        # If field is truly empty, ElementTree leaves off tag entirely which 
        # causes logbook parser to fail
        if text.text == "": text.text = " "
         
        xmlFile = open(fileName+'.xml',"w")
        rawString = ElementTree.tostring(log_entry, 'utf-8')
        
        # Adds newline after each closing tag
        parsedString = sub(r'(?=<[^/].*>)','\n',rawString) 
        
        xmlString = parsedString[1:]
        xmlFile.write(xmlString)
        
        # Close with newline so cron job parses correctly
        xmlFile.write("\n")  
        
        xmlFile.close()
        exporter = pg.exporters.ImageExporter(self.plot.plotItem)
        #exporter.parameters()['width'] = 550
        exporter.export(fileName+'.png')
        # PyQtGraph doesn't export PS files, so convert with linux 
        Popen('convert ' + fileName + '.png ' + fileName + '.ps', shell=True)
        
        sleep(0.3)
        Popen('convert ' + fileName + '.png -resize 500x500 ' + fileName
           + '.png', shell=True)
        sleep(0.35)
        path = "/u1/lcls/physics/logbook/data/"
        copy(fileName + '.ps', path)
        copy(fileName + '.png', path)
        copy(fileName + '.xml', path)
        self.statusBar().showMessage('Sent to LCLS Physics Logbook!', 10000)

    def MCCLog(self):
        exporter = pg.exporters.ImageExporter(self.plot.plotItem)
        exporter.export('/tmp/RTBSA.png')
        Popen("convert /tmp/RTBSA.png /tmp/RTBSA.ps", shell=True)
        sleep(0.1)
        printFile = "lpr -P" + 'elog_mcc' + " /tmp/RTBSA.ps"
        os.system(printFile)

    def stop(self):
        self.abort = True
        self.statusBar().showMessage('Stopped')
        self.ui.draw_button.setDisabled(False)
        QApplication.processEvents()
        
        try:
            self.val1pv.remove_callback(index = self.booyah)
            self.val1pv.disconnect()
        except:
            self.statusBar().showMessage('Stopped')
            
        try:
            self.val2pv.remove_callback(index = self.booyah2)
            self.val2pv.disconnect()
        except:
            self.statusBar().showMessage('Stopped')

    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        load_file_action = self.create_action("&Save plot",
            shortcut="Ctrl+S", slot=self.save_plot, 
            tip="Save the plot")
            
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (load_file_action, None, quit_action))
        
        self.help_menu = self.menuBar().addMenu("&Help")
        
        about_action = self.create_action("&About", shortcut = 'F1',
                  slot=self.on_about, tip = 'About')
                    
        self.add_actions(self.help_menu, (about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(self, text, slot=None, shortcut=None, icon=None, tip=None, 
             checkable=False, signal="triggered()"):
        action = QAction(text, self)
        
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action

    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.ui.widget.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)

    def on_about(self):
        msg = """Can you read this?  If so, congratulations. You are a magical, 
              marvelous troll."""
        QMessageBox.about(self, "About", msg.strip())

    def closeEvent(self,event):
        self.abort=True
        self.stop()

def main():
    app = QApplication(sys.argv)
    window = RTBSA()
    window.show()
    sys.exit(app.exec_())
    
if __name__ == "__main__":
    main()
