#!/usr/local/lcls/package/python/current/bin/python
# Written by Zimmer, modified by Lisa

import sys

from epics import caget, PV

import numpy as np
from numpy import polyfit, poly1d, polyval, corrcoef, std, mean
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from subprocess import CalledProcessError

from logbook import *

from rtbsa_ui import Ui_RTBSA as BSA_UI

from Constants import *


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

        # Generate list of BSA PVS
        try:
            updatedBSAPVs = check_output(['eget', '-ts', 'ds', '-a',
                                          'tag=LCLS.BSA.rootnames']).splitlines()[
                            1:-1]
            self.bsapvs = self.bsapvs + updatedBSAPVs

        # Backup for timeout error
        except CalledProcessError:
            self.bsapvs = self.bsapvs + bsapvs

        for pv in self.bsapvs:
            self.ui.listWidget.addItem(pv)
            self.ui.listWidget_2.addItem(pv)

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

        # Initial number of points
        self.numpoints = 2800

        # 20ms polling time
        self.updatetime = 20

        # Set initial polynomial fit to 2
        self.fitorder = 2

        self.dpi = 100

        self.ui.fitedit.setDisabled(True)
        self.ui.enter1.setDisabled(True)
        self.ui.enter2.setDisabled(True)
        self.ui.label.setDisabled(True)
        self.ui.listWidget.setDisabled(True)
        self.ui.listWidget_2.setDisabled(True)
        self.statusBar().showMessage('Hi there!  I missed you!')
        self.abort = True
        self.ui.parab_cb.setChecked(False)

        # Used to update plot
        self.timer = QTimer(self)

        self.rate = PV('EVNT:SYS0:1:LCLSBEAMRATE')
        self.menuBar().setStyleSheet(
            'QWidget{background-color:grey;color:purple}')
        self.create_menu()
        self.create_status_bar()

    def setUpGraph(self):
        self.plot = pg.PlotWidget(alpha=0.75)
        layout = QGridLayout()
        self.ui.widget.setLayout(layout)
        layout.addWidget(self.plot, 0, 0)
        self.plot.showGrid(1, 1)

    def loadStyleSheet(self):
        try:
            self.cssfile = "/home/physics/zimmerc/python/style.css"
            with open(self.cssfile, "r") as f:
                self.setStyleSheet(f.read())
        except:
            print "Error loading style sheet"
            pass

    def create_status_bar(self):
        self.status_text = QLabel()
        palette = QPalette()
        palette.setColor(palette.Foreground, Qt.magenta)
        self.statusBar().addWidget(self.status_text, 1)
        self.statusBar().setPalette(palette)

    def search(self, enter, widget):
        widget.clear()
        query = str(enter.text())
        for pv in self.bsapvs:
            if query.lower() in pv.lower():
                widget.addItem(pv)

    def search1(self):
        self.search(self.ui.enter1, self.ui.listWidget)

    def search2(self):
        self.search(self.ui.enter2, self.ui.listWidget_2)

    def setEnter(self, widget, enter, search, enter_rb):
        selection = widget.currentItem()
        enter.textChanged.disconnect()
        enter.setText(selection.text())
        QObject.connect(enter, SIGNAL("textChanged(const QString&)"), search)
        if not self.abort and enter_rb.isChecked():
            self.stop()
            self.on_draw()

    def setenter1(self):
        self.setEnter(self.ui.listWidget, self.ui.enter1, self.search1,
                 self.ui.enter1_rb)

    def setenter2(self):
        self.setEnter(self.ui.listWidget_2, self.ui.enter2, self.search2,
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
            self.pvlist.append(str(common.currentText() + 'HSTBR'))

        elif enter_rb.isChecked():
            if str(enter.text()).strip():
                self.pvlist.append(str(enter.text()) + 'HSTBR')
            else:
                self.statusBar().showMessage('Device ' + device
                                             + ' field blank. Aborting.', 10000)
                self.ui.draw_button.setEnabled(True)
                success = False

        return success

    def setValSynced(self):

        numBadShots = round((self.time2 - self.time1) * self.rate.value)

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

    def getData(self):
        # If i use Popen instead of pyepics caget, the program doesn't
        # start lagging if you change PVs (?!?!?). Some stupid bug in 
        # new pyepics.
        getdata = Popen("caget " + self.device, stdout=PIPE, shell=True)
        newdata = str(getdata.communicate()).split()[2:-1]
        newdata[-1] = newdata[-1][:-4]
        return [float(i) for i in newdata]

    def initialzeData(self):
        self.statusBar().showMessage('Initializing...')

        if self.ui.common1_rb.isChecked():
            self.device = str(self.ui.common1.currentText() + 'HSTBR')

        elif self.ui.enter1_rb.isChecked():
            if str(self.ui.enter1.text()).strip():
                self.device = str(self.ui.enter1.text() + 'HSTBR')
        else:
            return None

        return self.getData()

    def plotFit(self, newdata):
        self.curve = pg.PlotCurveItem(newdata[2800 - self.numpoints:2800],
                                      pen=1)
        self.plot.addItem(self.curve)
        self.plot.setTitle(self.device)
        self.xdata = range(self.numpoints)

        # Fit line
        if self.ui.line_cb.isChecked():
            self.getLinearFit(self.xdata, newdata[2800 - self.numpoints:2800],
                              False)
            self.plot.addItem(self.fit)

        # Fit polynomial
        elif self.ui.parab_cb.isChecked():
            self.ui.fitedit.setDisabled(False)
            self.getPolynomialFit(self.xdata,
                                  newdata[2800 - self.numpoints:2800],
                                  False)
            self.plot.addItem(self.parab)

    def genTimePlotA(self):
        newdata = self.initialzeData()

        # Using pyepics caget causes progressive slowdown of GUI!?!?

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
        self.timer = QTimer(self)
        self.timer.singleShot(self.updatetime, self.update_plot_HSTBR)

    def adjustVals(self):
        self.updateRate()
        numBadShots, val1synced, val2synced = self.setValSynced()

        blength = 2800 - numBadShots
        if (self.numpoints <= blength):
            self.val1 = val1synced[blength - self.numpoints:blength]
            self.val2 = val2synced[blength - self.numpoints:blength]
        else:
            self.val1 = val1synced
            self.val2 = val2synced

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

        self.booyah = self.val1pv.add_callback(self.ItDoneChanged)
        self.booyah2 = self.val2pv.add_callback(self.ItDoneChanged2)

        while (self.time1 == 1 or self.time2 == 2) and not self.abort:
            QApplication.processEvents()

        self.adjustVals()

    def getLinearFit(self, xdata, ydata, updateExistingPlot):
        try:
            (m, b) = polyfit(xdata, ydata, 1)
            fitdata = polyval([m, b], xdata)
            m = "{:.3e}".format(m)
            self.slope_text.setText('Slope: ' + str(m))
            if updateExistingPlot:
                self.fit.setData(xdata, fitdata)
            else:
                self.fit = pg.PlotCurveItem(xdata, fitdata, 'g-', linewidth=1)
        except:
            print "Error getting linear fit"
            pass

    def getPolynomialFit(self, xdata, ydata, updateExistingPlot):
        try:
            co = polyfit(xdata, ydata, self.fitorder)
            pol = poly1d(co)
            sorted1 = sorted(xdata)
            fit = pol(sorted1)

            if updateExistingPlot:
                self.parab.setData(sorted1, fit)
            else:
                self.parab = pg.PlotCurveItem(self.xdata, fit, pen=3)

            if self.fitorder == 2:
                self.slope_text.setText('Peak: ' + str(-co[1] / (2 * co[0])))

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

    # TODO This is near identical to plotFit - I'll ask Chris if the differences
    # are necessary
    def genABPlot(self):
        self.curve = pg.ScatterPlotItem(self.val1, self.val2, pen=1, symbol='x',
                                        size=5)
        self.plot.addItem(self.curve)
        self.plot.setTitle(self.pvlist[1] + ' vs. ' + self.pvlist[0])

        if self.ui.line_cb.isChecked():
            self.fit = pg.PlotCurveItem(pen=3)
            self.getLinearFit(self.val1, self.val2, True)
            self.plot.addItem(self.fit)

        elif self.ui.parab_cb.isChecked():
            self.parab = pg.PlotCurveItem(pen=3, size=2)
            self.getPolynomialFit(self.val1, self.val2, True)
            self.plot.addItem(self.parab)

    def genFFTPlot(self):
        newdata = self.initialzeData()
        newdata = newdata[2800 - self.numpoints:2800]
        newdata.extend(np.zeros(self.numpoints * 2).tolist())
        newdata = newdata - np.mean(newdata);
        ps = np.abs(np.fft.fft(newdata)) / len(newdata)
        self.FS = self.rate.value
        self.freqs = np.fft.fftfreq(len(newdata), 1.0 / self.FS)
        self.keep = self.freqs >= 0
        ps = ps[self.keep]
        self.freqs = self.freqs[self.keep]
        self.idx = np.argsort(self.freqs)
        self.curve = pg.PlotCurveItem(x=self.freqs[self.idx],
                                      y=ps[self.idx], pen=1)
        self.plot.addItem(self.curve)
        self.plot.setTitle(self.device)

    def genPlotAndSetTimer(self, genPlot):
        if self.abort:
            return False

        try:
            genPlot()
        except UnboundLocalError:
            self.statusBar().showMessage('No Data, Aborting Plotting Algorithm',
                                         10000)
            return False

        self.timer = QTimer(self)
        self.timer.singleShot(self.updatetime, self.update_BSA_Plot)
        self.statusBar().showMessage('Running')
        return True

    # Where the magic happens(well, where it starts to happen). This initializes
    # the BSA plotting and then starts a timer to update the plot.
    def on_draw(self):
        plotTypeIsValid = (self.ui.AvsT_cb.isChecked()
                           or self.ui.AvsB.isChecked()
                           or self.ui.AFFT.isChecked())

        if not plotTypeIsValid:
            self.statusBar().showMessage(
                'Pick a Plot Type (PV vs. time or B vs A)',
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

    def filterVals(self):
        # Filter out NaNs and ridiculous BLEN values
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

    def setPlotRanges(self):
        if self.ui.autoscale_cb.isChecked():
            mx = np.max(self.val2)
            mn = np.min(self.val2)

            if mn != mx:
                self.plot.setYRange(mn, mx)

            mx = np.max(self.val1)
            mn = np.min(self.val1)

            if mn != mx:
                self.plot.setXRange(mn, mx)

    def setPosAndText(self, attribute, value, posValX, posValY, textVal):
        value = "{:.3}".format(value)
        attribute.setPos(posValX, posValY)
        attribute.setText(textVal + str(value))

    def update_BSA_Plot(self):
        QApplication.processEvents()

        if self.abort:
            return

        self.plot.showGrid(self.ui.grid_cb.isChecked(),
                           self.ui.grid_cb.isChecked())

        self.adjustVals()
        self.filterVals()

        self.curve.setData(self.val1, self.val2)

        self.setPlotRanges()

        if self.ui.avg_cb.isChecked():
            a = mean(self.val2)
            self.setPosAndText(self.avg_text, a, min(self.val1),
                               min(self.val2), 'AVG: ')

        if self.ui.std_cb.isChecked():
            val1Min = min(self.val1)
            xPos = (val1Min + (val1Min + max(self.val1)) / 2) / 2
            s = std(self.val2)
            self.setPosAndText(self.std_text, s, xPos, min(self.val2), 'STD: ')

        if self.ui.corr_cb.isChecked():
            correlation = corrcoef(self.val1, self.val2)
            self.setPosAndText(self.corr_text, correlation, min(self.val1),
                               max(self.val2), "Corr. Coefficient: ")

        if self.ui.line_cb.isChecked():
            self.slope_text.setPos((min(self.val1) + max(self.val1)) / 2,
                                   min(self.val2))
            self.getLinearFit(self.val1, self.val2, True)

        elif self.ui.parab_cb.isChecked():
            self.slope_text.setPos((min(self.val1) + max(self.val1)) / 2,
                                   min(self.val2))
            self.getPolynomialFit(self.val1, self.val2, True)

        self.timer.singleShot(self.updatetime, self.update_BSA_Plot)

    def update_plot_HSTBR(self):

        self.plot.showGrid(self.ui.grid_cb.isChecked(),
                           self.ui.grid_cb.isChecked())

        QApplication.processEvents()

        if self.abort:
            return

        chopped = self.getData()[2800 - self.numpoints:2800]
        self.curve.setData(chopped)

        if self.ui.autoscale_cb.isChecked():
            mx = max(chopped)
            mn = min(chopped)
            if mx - mn > .00001:
                self.plot.setYRange(mn, mx)
                self.plot.setXRange(0, len(chopped))

        if self.ui.avg_cb.isChecked():
            self.setPosAndText(self.avg_text, mean(chopped), 0, min(chopped),
                               'AVG: ')

        if self.ui.std_cb.isChecked():
            self.setPosAndText(self.std_text, std(chopped), self.numpoints/4,
                               min(chopped), 'STD: ')

        if self.ui.corr_cb.isChecked():
            self.corr_text.setText('')

        if self.ui.line_cb.isChecked():
            self.slope_text.setPos(self.numpoints / 2, min(chopped))
            self.getLinearFit(self.xdata, chopped, True)

        elif self.ui.parab_cb.isChecked():
            self.slope_text.setPos(self.numpoints / 2, min(chopped))
            self.getPolynomialFit(self.xdata, chopped, True)

        self.timer.singleShot(40, self.update_plot_HSTBR)

    def update_plot_FFT(self):
        self.plot.showGrid(self.ui.grid_cb.isChecked(),
                           self.ui.grid_cb.isChecked())
        QApplication.processEvents()

        if self.abort:
            return

        newdata = self.getData()
        newdata = newdata[2800 - self.numpoints:2800]
        newdata = np.array(newdata)
        nans, x = np.isnan(newdata), lambda z: z.nonzero()[0]
        # interpolate nans
        newdata[nans] = np.interp(x(nans), x(~nans), newdata[~nans])
        newdata = newdata - np.mean(newdata);
        newdata = newdata.tolist()
        newdata.extend(np.zeros(self.numpoints * 2).tolist())
        ps = np.abs(np.fft.fft(newdata)) / len(newdata)
        self.FS = self.rate.value
        self.freqs = np.fft.fftfreq(len(newdata), 1.0 / self.FS)
        self.keep = (self.freqs >= 0)
        ps = ps[self.keep]
        self.freqs = self.freqs[self.keep]
        self.idx = np.argsort(self.freqs)
        self.curve.setData(x=self.freqs[self.idx], y=ps[self.idx])

        if self.ui.autoscale_cb.isChecked():
            mx = max(ps)
            mn = min(ps)
            if mx - mn > .00001:
                self.plot.setYRange(mn, mx)
                self.plot.setXRange(min(self.freqs), max(self.freqs))

        self.timer.singleShot(40, self.update_plot_FFT)

    # Callback function for PV1
    def ItDoneChanged(self, pvname=None, value=None, timestamp=None, **kw):
        self.time1 = timestamp
        self.val1pre = value

    # Callback function for PV2
    def ItDoneChanged2(self, pvname=None, value=None, timestamp=None, **kw):
        self.time2 = timestamp
        self.val2pre = value

    def AvsTClick(self):
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
            self.timer.singleShot(150, self.on_draw)

            # Mechanism to fix a race condition- start a timer which allows the 
            # main loop to regain control and see that the abort variable has 
            # turned positive, and then the timer restarts the drawing process 
            # after a short time, giving the active loop time to stop.  This 
            # fixes a problem I was seeing (Gibbs noticed first) with a double 
            # plot showing up (rare occurence, only happened if the user tried 
            # to plot a blank PV name from the enter line)

            # QObject.connect(self.timer2, SIGNAL("timeout()"), self.on_draw)
            # self.timer2.start(self.updatetime)

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
            self.statusBar().showMessage('Enter an integer, 1-10', 6000)
            return

        if self.fitorder > 10 or self.fitorder < 1:
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
            # Setup for single PV plotting
            if self.ui.AvsT_cb.isChecked():
                self.genTimePlotA()

            elif self.ui.AvsB.isChecked():
                self.genABPlot()
            else:
                newdata = self.getData()
                newdata = newdata[2800 - self.numpoints:2800]
                newdata.extend(np.zeros(self.numpoints * 2).tolist())
                newdata = newdata - np.mean(newdata);
                ps = np.abs(np.fft.fft(newdata)) / len(newdata)
                rate = {4: 1, 5: 10, 6: 30, 7: 60, 8: 120}
                i = caget('IOC:BSY0:MP01:BYKIK_RATE')
                self.FS = rate[i];
                self.freqs = np.fft.fftfreq(self.numpoints, 1.0 / self.FS)
                self.keep = (self.freqs >= 0)
                self.freqs = self.freqs[self.keep]
                ps = ps[self.keep]
                self.idx = np.argsort(self.freqs)
                self.curve = pg.PlotCurveItem(x=self.freqs[self.idx],
                                              y=ps[self.idx], pen=1)
                self.plot.addItem(self.curve)
                self.plot.setTitle(self.device)

        except:
            print "Error reinitializing plot"
            pass

    def logbook(self):
        logbook('Python Real-Time BSA', 'BSA Data',
                str(self.numpoints) + ' points', self.plot.plotItem)
        self.statusBar().showMessage('Sent to LCLS Physics Logbook!', 10000)

    def MCCLog(self):
        MCCLog('/tmp/RTBSA.png', '/tmp/RTBSA.ps', self.plot.plotItem)

    def stop(self):
        self.abort = True
        self.statusBar().showMessage('Stopped')
        self.ui.draw_button.setDisabled(False)
        QApplication.processEvents()

        try:
            self.val1pv.remove_callback(index=self.booyah)
            self.val1pv.disconnect()
        except:
            self.statusBar().showMessage('Stopped')

        try:
            self.val2pv.remove_callback(index=self.booyah2)
            self.val2pv.disconnect()
        except:
            self.statusBar().showMessage('Stopped')

    def create_menu(self):
        self.file_menu = self.menuBar().addMenu("&File")

        load_file_action = self.create_action("&Save plot",
                                              shortcut="Ctrl+S",
                                              slot=self.save_plot,
                                              tip="Save the plot")

        quit_action = self.create_action("&Quit", slot=self.close,
                                         shortcut="Ctrl+Q",
                                         tip="Close the application")

        self.add_actions(self.file_menu,
                         (load_file_action, None, quit_action))

        self.help_menu = self.menuBar().addMenu("&Help")

        about_action = self.create_action("&About", shortcut='F1',
                                          slot=self.on_about, tip='About')

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


def main():
    app = QApplication(sys.argv)
    window = RTBSA()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()