#!/usr/local/lcls/package/python/current/bin/python
# Written by Zimmer, refactored by Lisa

import sys
import time

from epics import PV

import numpy as np
from numpy import polyfit, poly1d, polyval, corrcoef, std, mean
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from subprocess import CalledProcessError, check_output, PIPE
from itertools import compress

from logbook import *
from rtbsa_ui import Ui_RTBSA
from Constants import *

class RTBSA(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.ui = Ui_RTBSA()
        self.ui.setupUi(self)
        self.setWindowTitle('Real Time BSA')
        self.loadStyleSheet()
        self.setUpGraph()

        # enter 1 is the text input box for device A, and 2 is for B
        QObject.connect(self.ui.enter1, SIGNAL("textChanged(const QString&)"),
                        self.searchA)
        QObject.connect(self.ui.enter2, SIGNAL("textChanged(const QString&)"),
                        self.searchB)

        # Changes the text in the input box to match the selection from the list
        self.ui.listWidget.itemClicked.connect(self.setEnterA)
        self.ui.listWidget_2.itemClicked.connect(self.setEnterB)

        self.bsapvs = ['GDET:FEE1:241:ENRC', 'GDET:FEE1:242:ENRC',
                       'GDET:FEE1:361:ENRC', 'GDET:FEE1:362:ENRC']

        # Generate list of BSA PVS
        try:
            BSAPVs = check_output(['eget', '-ts', 'ds', '-a',
                                   'tag=LCLS.BSA.rootnames']).splitlines()[1:-1]
            self.bsapvs = self.bsapvs + BSAPVs

        # Backup for eget error
        except CalledProcessError:
            print "Unable to pull most recent PV list"
            self.bsapvs = self.bsapvs + bsapvs

        for pv in self.bsapvs:
            self.ui.listWidget.addItem(pv)
            self.ui.listWidget_2.addItem(pv)

        # Dropdown menu for device A (add common BSA PV's and make bunch length
        # the default selection)
        self.ui.common1.addItems(commonlist)
        self.ui.common1.setCurrentIndex(21)
        self.ui.common1.activated.connect(self.inputActivated)

        # Dropdown menu for device B
        self.ui.common2.addItems(commonlist)
        self.ui.common2.activated.connect(self.inputActivated)

        # All the checkboxes in the Settings section
        self.ui.AvsT_cb.clicked.connect(self.AvsTClick)
        self.ui.AvsB.clicked.connect(self.AvsBClick)
        self.ui.AFFT.clicked.connect(self.AFFTClick)
        self.ui.avg_cb.clicked.connect(self.avg_click)
        self.ui.std_cb.clicked.connect(self.std_click)
        self.ui.corr_cb.clicked.connect(self.corr_click)
        self.ui.parab_cb.clicked.connect(self.parab_click)
        self.ui.line_cb.clicked.connect(self.line_click)

        # All the buttons in the Controls section
        self.ui.draw_button.clicked.connect(self.on_draw)
        self.ui.stop_button.clicked.connect(self.stop)
        self.ui.log_button.clicked.connect(self.logbook)
        self.ui.mcclog_button.clicked.connect(self.MCCLog)

        # fitedit is the text input box for "Order"
        self.ui.fitedit.returnPressed.connect(self.fitOrderActivated)

        # The radio buttons that enable the dropdown menus
        self.ui.common1_rb.clicked.connect(self.common_1_click)
        self.ui.common2_rb.clicked.connect(self.common_2_click)

        # The radio buttons that enable the search bars
        self.ui.enter1_rb.clicked.connect(self.enter_1_click)
        self.ui.enter2_rb.clicked.connect(self.enter_2_click)

        # TODO ask Chris if we can remove this functionality
        # Triggers a redrawing upon pressing enter in the search bar.
        # Proper usage should be using the search bar to search, and selecting
        # from the results in the list. If it's not in the list, it's an invalid
        # PV with no reason to attempt plotting
        # self.ui.enter1.returnPressed.connect(self.inputActivated)
        # self.ui.enter2.returnPressed.connect(self.inputActivated)

        # Pressing enter in the text input boxes for points and std dev triggers
        # updating the plot
        self.ui.points.returnPressed.connect(self.points_entered)
        self.ui.numStdDevs.returnPressed.connect(self.stdDevEntered)

        # Initial number of points
        self.numpoints = 2800

        # Initial number of standard deviations
        self.stdDevstoKeep = 3.0

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
        self.menuBar().setStyleSheet('QWidget{background-color:grey;color:purple}')
        self.create_menu()
        self.create_status_bar()

        # A dictionary to store the PV names, initialized to None
        self.devices = {"A": None, "B": None}

        # Versions of data buffers A and B that are filtered by standard
        # deviation. Didn't want to edit those buffers directly so that we could
        # unfliter or refilter with a different number more efficiently
        self.filteredBufferA = None
        self.filteredBufferB = None

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

    # Effectively an autocomplete
    def search(self, enter, widget):
        widget.clear()
        query = str(enter.text())
        for pv in self.bsapvs:
            if query.lower() in pv.lower():
                widget.addItem(pv)

    def searchA(self):
        self.search(self.ui.enter1, self.ui.listWidget)

    def searchB(self):
        self.search(self.ui.enter2, self.ui.listWidget_2)

    def setEnter(self, widget, enter, search, enter_rb):
        selection = widget.currentItem()
        enter.textChanged.disconnect()
        enter.setText(selection.text())
        QObject.connect(enter, SIGNAL("textChanged(const QString&)"), search)
        if not self.abort and enter_rb.isChecked():
            self.stop()
            self.on_draw()

    def setEnterA(self):
        self.setEnter(self.ui.listWidget, self.ui.enter1, self.searchA,
                      self.ui.enter1_rb)

    def setEnterB(self):
        self.setEnter(self.ui.listWidget_2, self.ui.enter2, self.searchB,
                      self.ui.enter2_rb)

    def correctInput(self, errorMessage, acceptableTxt, textBox):
        self.statusBar().showMessage(errorMessage, 6000)
        textBox.setText(acceptableTxt)

    def correctNumpoints(self, errorMessage, acceptableValue):
        self.correctInput(errorMessage, str(acceptableValue), self.ui.points)
        self.numpoints = acceptableValue

    def correctStdDevs(self, errorMessage, acceptableValue):
        self.correctInput(errorMessage, str(acceptableValue), self.ui.numStdDevs)
        self.stdDevstoKeep = acceptableValue

    def stdDevEntered(self):
        try:
            self.stdDevstoKeep = float(self.ui.numStdDevs.text())
        except ValueError:
            self.correctStdDevs('Enter a float > 0', 3.0)
            return

        # Is there a way to combine an except and an if?
        if self.stdDevstoKeep <= 0:
            self.correctStdDevs('Enter a float > 0', 3.0)
            return

    def points_entered(self):
        try:
            self.numpoints = int(self.ui.points.text())
        except ValueError:
            self.correctNumpoints('Enter an integer, 1 to 2800', 120)
            return

        if self.numpoints > 2800:
            self.correctNumpoints('Max # points is 2800', 2800)
            return

        if self.numpoints < 1:
            self.correctNumpoints('Min # points is 1', 1)
            return

        self.reinitialize_plot()

    def populateDevices(self, common_rb, common, enter_rb, enter, device):

        # HSTBR is what gets the beam synchronous data buffer
        if common_rb.isChecked():
            self.devices[device] = str(common.currentText()) + 'HSTBR'

        elif enter_rb.isChecked():
            pv = str(enter.text()).strip()
            # Checks that it's non empty and that it's a BSA pv
            if pv and pv in self.bsapvs:
                self.devices[device] = pv + 'HSTBR'
            else:
                self.statusBar().showMessage('Device ' + device
                                             + ' invalid. Aborting.', 10000)
                self.ui.draw_button.setEnabled(True)
                return False

        return True

    ############################################################################
    # Time 1 is when Device A started acquiring data, and Time 2 is when Device
    # B started acquiring data. Since they're not guaranteed to start
    # acquisition at the same time, one data buffer might be ahead of the other,
    # meaning that the intersection of the two buffers would not include the
    # first n elements of one and the last n elements of the other. See the
    # diagram below, where the dotted line represents the time axis (one buffer
    # is contained  by square brackets [], the other by curly braces {}, and the
    # times where each starts  and ends is indicated right underneath).
    #
    #
    #          [           {                            ]           }
    # <----------------------------------------------------------------------> t
    #       t1_start    t2_start                     t1_end      t2_end
    #
    #
    # Note that both buffers are of the same size (self.numpoints) so that:
    # (t1_end - t1_start) = (t2_end - t2_start)
    #
    # From the diagram, we see that only the time between t2_start and t1_end
    # contains data from both buffers (t1_start to t2_start only contains data
    # from buffer 1, and t1_end to t2_end only contains data from buffer 2).
    # Using that, we can chop the beginning of buffer 1 and the end of buffer 2
    # so that we're only left with the overlapping region.
    #
    # In order to figure out how many points we need to chop from each buffer
    # (it's the same number for both since they're both the same size), we
    # multiply the time delta by the beam rate (yay dimensional analysis!):
    # seconds * (shots/second) = (number of shots)
    ############################################################################
    def setValSynced(self):

        numBadShots = round((self.startTimeB - self.startTimeA) * self.rate.value)

        # Gets opposite indices depending on which time is greater (and [0:2800]
        # if they're equal)
        self.dataBufferA = self.unsynchedBufferA[max(0, numBadShots)
                                                 :min(2800, 2800 + numBadShots)]
        self.dataBufferB = self.unsynchedBufferB[max(0, -numBadShots)
                                                 :min(2800, 2800 - numBadShots)]

        return abs(numBadShots)

    # A spin loop that accounts for intermediate values when changing rates
    def updateRate(self):
        # self.rate is a PV, such that .value is shorthand for .getval
        rate = self.rate.value

        start_time = time.time()
        stuck = False
        while rate not in [120.0, 60.0, 30.0, 10.0]:
            QApplication.processEvents()
            rate = self.rate.value
            if time.time() - start_time > 1:
                stuck = True
                self.statusBar().showMessage("Waiting for beam rate to be 10, "
                                             "30, 60, or 120 Hz...")

        if stuck:
            self.statusBar().showMessage("Beam rate at allowed value")
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
        getdata = Popen("caget " + self.fullPVName, stdout=PIPE, shell=True)

        # TODO What are these indices?
        newdata = str(getdata.communicate()).split()[2:-1]
        newdata[-1] = newdata[-1][:-4]
        return [float(i) for i in newdata]

    def initializeData(self):
        self.statusBar().showMessage('Initializing...')

        if self.ui.common1_rb.isChecked():
            self.fullPVName = str(self.ui.common1.currentText() + 'HSTBR')

        elif self.ui.enter1_rb.isChecked():
            pv = str(self.ui.enter1.text()).strip()
            if pv and pv in self.bsapvs:
                self.fullPVName = pv + 'HSTBR'
            else:
                return None
        else:
            return None

        return self.getData()

    def adjustVals(self):
        self.updateRate()
        numBadShots = self.setValSynced()

        blength = 2800 - numBadShots

        # Make sure the buffer size doesn't exceed the desired number of points
        if (self.numpoints < blength):
            self.dataBufferA = self.dataBufferA[blength - self.numpoints:blength]
            self.dataBufferB = self.dataBufferB[blength - self.numpoints:blength]

    def updateValsFromInput(self):

        if not self.populateDevices(self.ui.common1_rb, self.ui.common1,
                                    self.ui.enter1_rb, self.ui.enter1, "A"):
            return False

        if not self.populateDevices(self.ui.common2_rb, self.ui.common2,
                                    self.ui.enter2_rb, self.ui.enter2, "B"):
            return False

        self.dataBufferA, self.dataBufferB = [], []

        # Without the time parameter, we wouldn't get the timestamp
        self.val1pv = PV(self.devices["A"], form='time')
        self.val2pv = PV(self.devices["B"], form='time')

        self.statusBar().showMessage('Initializing/Syncing (be patient, '
                                     + 'may take 5 seconds)...')

        self.startTimeA, self.startTimeB = None, None

        # For some reason, callbacks are the only way to get a timestamp
        self.val1CallbackIndex = self.val1pv.add_callback(self.val1Callback)
        self.val2CallbackIndex = self.val2pv.add_callback(self.val2Callback)

        while (not self.startTimeA or not self.startTimeB) and not self.abort:
            QApplication.processEvents()

        self.adjustVals()
        return True

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
            xdataSorted = sorted(xdata)
            fit = pol(xdataSorted)

            if updateExistingPlot:
                self.parab.setData(xdataSorted, fit)
            else:
                self.parab = pg.PlotCurveItem(xdataSorted, fit, pen=3, size=2)

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

    def plotFit(self, xdata, ydata, title):
        self.plot.addItem(self.curve)
        self.plot.setTitle(title)

        # Fit line
        if self.ui.line_cb.isChecked():
            self.getLinearFit(xdata, ydata, False)
            self.plot.addItem(self.fit)

        # Fit polynomial
        elif self.ui.parab_cb.isChecked():
            self.ui.fitedit.setDisabled(False)
            self.getPolynomialFit(xdata, ydata, False)
            self.plot.addItem(self.parab)

    def genTimePlotA(self):
        newdata = self.initializeData()

        if not newdata:
            self.statusBar().showMessage('Invalid PV? Unable to get data.'
                                         + ' Aborting.', 10000)
            self.ui.draw_button.setEnabled(True)
            return

        self.curve = pg.PlotCurveItem(newdata[2800 - self.numpoints:2800],
                                      pen=1)
        self.plot.addItem(self.curve)

        self.xdata = range(self.numpoints)

        self.plotFit(self.xdata, newdata[2800 - self.numpoints:2800],
                     self.fullPVName)

    def plotCurveAndFit(self, xdata, ydata):
        self.curve = pg.ScatterPlotItem(xdata, ydata, pen=1, symbol='x', size=5)
        self.plot.addItem(self.curve)
        self.plotFit(xdata, ydata,
                     self.devices["B"] + ' vs. ' + self.devices["A"])

    def genABPlot(self):
        if self.ui.filterByStdDevs.isChecked():
            self.plotCurveAndFit(self.filteredBufferA, self.filteredBufferB)
        else:
            self.plotCurveAndFit(self.dataBufferA, self.dataBufferB)

    # TODO I have no idea what's happening here
    def updateFFTPlot(self, newdata, updateExistingPlot):

        if not newdata:
            return None

        newdata = newdata[2800 - self.numpoints:2800]

        newdata = np.array(newdata)
        nans, x = np.isnan(newdata), lambda z: z.nonzero()[0]
        # interpolate nans
        newdata[nans] = np.interp(x(nans), x(~nans), newdata[~nans])
        # remove DC component
        newdata = newdata - np.mean(newdata)
        newdata = newdata.tolist()

        newdata.extend(np.zeros(self.numpoints * 2).tolist())

        ps = np.abs(np.fft.fft(newdata)) / len(newdata)

        self.FS = self.updateRate()

        self.freqs = np.fft.fftfreq(len(newdata), 1.0 / self.FS)
        self.keep = (self.freqs >= 0)
        ps = ps[self.keep]
        self.freqs = self.freqs[self.keep]
        self.idx = np.argsort(self.freqs)

        if updateExistingPlot:
            self.curve.setData(x=self.freqs[self.idx], y=ps[self.idx])
        else:
            self.curve = pg.PlotCurveItem(x=self.freqs[self.idx],
                                          y=ps[self.idx], pen=1)
        self.plot.addItem(self.curve)
        self.plot.setTitle(self.fullPVName)

        return ps


    def genFFTPlot(self):
        self.updateFFTPlot(self.initializeData(), False)

    def genPlotAndSetTimer(self, genPlot, updateMethod):
        if self.abort:
            return

        try:
            genPlot()
        except UnboundLocalError:
            self.statusBar().showMessage('No Data, Aborting Plotting Algorithm',
                                         10000)
            return

        self.timer = QTimer(self)
        self.timer.singleShot(self.updatetime, updateMethod)
        self.statusBar().showMessage('Running')

    # Where the magic happens (well, where it starts to happen). This
    # initializes the BSA plotting and then starts a timer to update the plot.
    def on_draw(self):
        plotTypeIsValid = (self.ui.AvsT_cb.isChecked()
                           or self.ui.AvsB.isChecked()
                           or self.ui.AFFT.isChecked())

        if not plotTypeIsValid:
            self.statusBar().showMessage('Pick a Plot Type (PV vs. time '
                                         'or B vs A)', 10000)
            return

        self.ui.draw_button.setDisabled(True)
        self.abort = False

        self.cleanPlot()

        ####Plot history buffer for one PV####
        if self.ui.AvsT_cb.isChecked():
            if self.populateDevices(self.ui.common1_rb, self.ui.common1,
                                    self.ui.enter1_rb, self.ui.enter1, "A"):
                self.genPlotAndSetTimer(self.genTimePlotA,
                                        self.update_plot_HSTBR)

        ####Plot for 2 PVs####
        elif self.ui.AvsB.isChecked():
            if self.updateValsFromInput():
                self.genPlotAndSetTimer(self.genABPlot, self.update_BSA_Plot)

        ####Plot power spectrum####
        else:
            if self.populateDevices(self.ui.common1_rb, self.ui.common1,
                                    self.ui.enter1_rb, self.ui.enter1, "A"):
                self.genPlotAndSetTimer(self.genFFTPlot, self.update_plot_FFT)

    # Need to filter out errant indices from both buffers to keep them
    # synchronized
    def filterData(self, dataBuffer, filterFunc, changeOriginal):
        mask = [filterFunc(x) for x in dataBuffer]

        if changeOriginal:
            self.dataBufferA = list(compress(self.dataBufferA, mask))
            self.dataBufferB = list(compress(self.dataBufferB, mask))
        else:
            self.filteredBufferA = list(compress(self.dataBufferA, mask))
            self.filteredBufferB = list(compress(self.dataBufferB, mask))

    # This PV gets insane values, apparently
    def filterPeakCurrent(self):
        filterFunc = lambda x: x < 12000
        if self.devices["A"] == "BLEN:LI24:886:BIMAXHSTBR":
            self.filterData(self.dataBufferA, filterFunc, True)
        if self.devices["B"] == "BLEN:LI24:886:BIMAXHSTBR":
            self.filterData(self.dataBufferB, filterFunc, True)

    def filterNans(self):
        filterFunc = lambda x: not np.isnan(x)
        self.filterData(self.dataBufferA, filterFunc, True)
        self.filterData(self.dataBufferB, filterFunc, True)

    def StdDevFilterFunc(self, average, stdDev):
        return lambda x: abs(x - average) < self.stdDevstoKeep*stdDev

    def filterStdDev(self):

        self.filterData(self.dataBufferA,
                        self.StdDevFilterFunc(mean(self.dataBufferA),
                                              std(self.dataBufferA)),
                        False)
        self.filterData(self.dataBufferB,
                        self.StdDevFilterFunc(mean(self.dataBufferB),
                                              std(self.dataBufferB)),
                        False)

    def setPlotRanges(self, bufferA, bufferB):
        mx = np.max(bufferB)
        mn = np.min(bufferB)

        if mn != mx:
            self.plot.setYRange(mn, mx)

        mx = np.max(bufferA)
        mn = np.min(bufferA)

        if mn != mx:
            self.plot.setXRange(mn, mx)

    def setPosAndText(self, attribute, value, posValX, posValY, textVal):
        value = "{:.3}".format(value)
        attribute.setPos(posValX, posValY)
        attribute.setText(textVal + str(value))

    def updateLabelsAndFit(self, bufferA, bufferB):
        self.curve.setData(bufferA, bufferB)

        self.setPlotRanges(bufferA, bufferB)

        # Logic to determine positions of labels when not running autoscale
        if self.ui.avg_cb.isChecked():
            self.setPosAndText(self.avg_text, mean(bufferB), min(bufferA),
                               min(bufferB), 'AVG: ')

        if self.ui.std_cb.isChecked():
            val1Min = min(bufferA)
            xPos = (val1Min + (val1Min + max(bufferA)) / 2) / 2

            self.setPosAndText(self.std_text, std(bufferB), xPos, min(bufferB),
                               'STD: ')

        if self.ui.corr_cb.isChecked():
            correlation = corrcoef(bufferA, bufferB)
            self.setPosAndText(self.corr_text, correlation, min(bufferA),
                               max(bufferB), "Corr. Coefficient: ")

        if self.ui.line_cb.isChecked():
            self.slope_text.setPos((min(bufferA) + max(bufferA)) / 2,
                                   min(bufferB))
            self.getLinearFit(bufferA, bufferB, True)

        elif self.ui.parab_cb.isChecked():
            self.slope_text.setPos((min(bufferA) + max(bufferA)) / 2,
                                   min(bufferB))
            self.getPolynomialFit(bufferA, bufferB, True)

    def update_BSA_Plot(self):
        QApplication.processEvents()

        if self.abort:
            return

        self.plot.showGrid(self.ui.grid_cb.isChecked(),
                           self.ui.grid_cb.isChecked())

        self.adjustVals()
        self.filterNans()
        self.filterPeakCurrent()

        if self.ui.filterByStdDevs.isChecked():
            self.filterStdDev()
            self.updateLabelsAndFit(self.filteredBufferA, self.filteredBufferB)
        else:
            self.updateLabelsAndFit(self.dataBufferA, self.dataBufferB)

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

        ps = self.updateFFTPlot(self.getData(), True)

        if self.ui.autoscale_cb.isChecked():
            mx = max(ps)
            mn = min(ps)
            if mx - mn > .00001:
                self.plot.setYRange(mn, mx)
                self.plot.setXRange(min(self.freqs), max(self.freqs))

        self.timer.singleShot(40, self.update_plot_FFT)

    # Callback function for Device A
    # value is the buffer because we're monitoring the HSTBR PV
    def val1Callback(self, pvname=None, value=None, timestamp=None, **kw):
        self.startTimeA = timestamp
        self.unsynchedBufferA = value

    # Callback function for Device B
    def val2Callback(self, pvname=None, value=None, timestamp=None, **kw):
        self.startTimeB = timestamp
        self.unsynchedBufferB = value

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
        self.inputActivated()

    def inputActivated(self):
        if not self.abort:
            self.stop()
            self.timer.singleShot(150, self.on_draw)

    def common_2_click(self):
        if self.ui.common2_rb.isChecked():
            self.ui.common2.setEnabled(True)
            self.ui.enter2_rb.setChecked(False)
            self.ui.enter2.setDisabled(True)
            self.ui.listWidget_2.setDisabled(True)
        else:
            self.ui.common2.setEnabled(False)
        self.inputActivated()

    def line_click(self):
        self.ui.parab_cb.setChecked(False)
        self.ui.fitedit.setDisabled(True)
        self.ui.label.setDisabled(True)
        self.reinitialize_plot()

    def fitOrderActivated(self):
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
                self.updateFFTPlot(self.getData(), False)

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
            self.val1pv.remove_callback(index=self.val1CallbackIndex)
            self.val1pv.disconnect()
        except:
            self.statusBar().showMessage('Stopped')

        try:
            self.val2pv.remove_callback(index=self.val2CallbackIndex)
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
        msg = ("Can you read this?  If so, congratulations. You are a magical, "
              + "marvelous troll.")
        QMessageBox.about(self, "About", msg.strip())


def main():
    app = QApplication(sys.argv)
    window = RTBSA()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
