#!/usr/local/lcls/package/python/current/bin/python
# Written by Zimmer, edited by Ahmed, refactored by Lisa
from os import path
from sys import argv, exit
from time import time

from epics import PV

# TODO import these with the namespace
from numpy import (polyfit, poly1d, polyval, corrcoef, std, mean, concatenate,
                   empty, nan, zeros, isnan, linalg, abs,
                   fft, argsort, interp, arange, nanmin, nanmax)

from PyQt4.QtCore import QTimer, QObject, SIGNAL, Qt
from PyQt4.QtGui import (QMainWindow, QLabel, QGridLayout, QPalette,
                         QApplication, QAction, QFileDialog, QIcon, QMessageBox)
from pyqtgraph import PlotWidget, PlotCurveItem, ScatterPlotItem, TextItem

from scipy.stats import nanmean, nanstd
from subprocess import CalledProcessError, check_output

from rtbsa_UI import Ui_RTBSA
import rtbsaUtils


# noinspection PyArgumentList,PyCompatibility
class RTBSA(QMainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.help_menu = self.menuBar().addMenu("&Help")
        self.file_menu = self.menuBar().addMenu("&File")
        self.status_text = QLabel()
        self.plot = PlotWidget(alpha=0.75)
        self.ui = Ui_RTBSA()
        self.ui.setupUi(self)
        self.setWindowTitle('Real Time BSA')
        self.loadStyleSheet()
        self.setUpGraph()

        self.bsapvs = ['GDET:FEE1:241:ENRC', 'GDET:FEE1:242:ENRC',
                       'GDET:FEE1:361:ENRC', 'GDET:FEE1:362:ENRC']

        self.populateBSAPVs()
        self.connectGuiFunctions()

        # Initial number of points
        self.numPoints = 2800

        # Initial number of standard deviations
        self.stdDevstoKeep = 3.0

        # 20ms polling time
        self.updateTime = 50

        # Set initial polynomial fit to 2
        self.fitOrder = 2

        self.disableInputs()
        self.abort = True

        # Used to update plot
        self.timer = QTimer(self)

        self.ratePV = PV('IOC:IN20:EV01:RG01_ACTRATE')

        self.menuBar().setStyleSheet('QWidget{background-color:grey;color:purple}')
        self.create_menu()
        self.create_status_bar()

        # The PV names
        self.devices = {"A": "", "B": ""}

        self.pvObjects = {"A": None, "B": None}

        # The raw, unsynchronized, unfiltered buffers
        self.rawBuffers = {"A": empty(2800), "B": empty(2800)}

        # The times when each buffer finished its last data acquisition
        self.timeStamps = {"A": None, "B": None}

        self.synchronizedBuffers = {"A": empty(2800), "B": empty(2800)}

        # Versions of data buffers A and B that are filtered by standard
        # deviation. Didn't want to edit those buffers directly so that we could
        # unfilter or refilter with a different number more efficiently
        self.filteredBuffers = {"A": empty(2800), "B": empty(2800)}

        # Text objects that appear on the plot
        self.text = {"avg": None, "std": None, "slope": None, "corr": None}

        # All things plot related!
        self.plotAttributes = {"curve": None, "fit": None, "parab": None,
                               "frequencies": None}

        # Used for the kill swtich
        self.counter = {"A": 0, "B": 0}

        # Used to implement scrolling for time plots
        self.currIdx = {"A": 0, "B": 0}

    def getRate(self):
        return rtbsaUtils.rateDict[self.ratePV.value]

    def disableInputs(self):
        self.ui.fitOrder.setDisabled(True)
        self.ui.searchInputA.setDisabled(True)
        self.ui.searchInputB.setDisabled(True)
        self.ui.labelFitOrder.setDisabled(True)
        self.ui.bsaListA.setDisabled(True)
        self.ui.bsaListB.setDisabled(True)
        self.statusBar().showMessage('Hi there!  I missed you!')
        self.ui.checkBoxPolyFit.setChecked(False)

    def populateBSAPVs(self):
        # Generate list of BSA PVS
        try:
            BSAPVs = check_output(['eget', '-ts', 'ds', '-a',
                                   'tag=LCLS.BSA.rootnames']).splitlines()[1:-1]
            self.bsapvs.extend(BSAPVs)

        # Backup for eget error
        except CalledProcessError:
            print("Unable to pull most recent PV list")
            # bsaPVs is pulled from the Constants file
            self.bsapvs.extend(rtbsaUtils.bsaPVs)

        for pv in self.bsapvs:
            self.ui.bsaListA.addItem(pv)
            self.ui.bsaListB.addItem(pv)

    def connectGuiFunctions(self):
        # enter 1 is the text input box for device A, and 2 is for B
        QObject.connect(self.ui.searchInputA, SIGNAL("textChanged(const QString&)"),
                        self.searchA)
        QObject.connect(self.ui.searchInputB, SIGNAL("textChanged(const QString&)"),
                        self.searchB)

        # Changes the text in the input box to match the selection from the list
        self.ui.bsaListA.itemClicked.connect(self.setEnterA)
        self.ui.bsaListB.itemClicked.connect(self.setEnterB)

        # Dropdown menu for device A (add common BSA PV's)
        self.ui.dropdownA.addItems(rtbsaUtils.commonlist)

        # Make bunch length default selection
        index = rtbsaUtils.commonlist.index("BLEN:LI24:886:BIMAX")
        self.ui.dropdownA.setCurrentIndex(index)

        self.ui.dropdownA.activated.connect(self.inputActivated)

        # Dropdown menu for device B
        self.ui.dropdownB.addItems(rtbsaUtils.commonlist)
        self.ui.dropdownB.activated.connect(self.inputActivated)

        # All the checkboxes in the Settings section
        self.ui.checkBoxAvsT.clicked.connect(self.AvsTClick)
        self.ui.checkBoxBvsA.clicked.connect(self.AvsBClick)
        self.ui.checkBoxFFT.clicked.connect(self.AFFTClick)
        self.ui.checkBoxShowAve.clicked.connect(self.avg_click)
        self.ui.checkBoxShowStdDev.clicked.connect(self.std_click)
        self.ui.checkBoxCorrCoeff.clicked.connect(self.corr_click)
        self.ui.checkBoxPolyFit.clicked.connect(self.parab_click)
        self.ui.checkBoxLinFit.clicked.connect(self.line_click)
        self.ui.checkBoxShowGrid.clicked.connect(self.showGrid)

        # All the buttons in the Controls section
        self.ui.startButton.clicked.connect(self.initializePlot)
        self.ui.stopButton.clicked.connect(self.stop)
        self.ui.lclsLogButton.clicked.connect(self.logbook)
        self.ui.mccLogButton.clicked.connect(self.MCCLog)

        # fitedit is the text input box for "Order"
        self.ui.fitOrder.returnPressed.connect(self.fitOrderActivated)

        # The radio buttons that enable the dropdown menus
        self.ui.dropdownButtonA.clicked.connect(self.common_1_click)
        self.ui.dropdownButtonB.clicked.connect(self.common_2_click)

        # The radio buttons that enable the search bars
        self.ui.searchButtonA.clicked.connect(self.enter_1_click)
        self.ui.searchButtonB.clicked.connect(self.enter_2_click)

        # Pressing enter in the text input boxes for points and std dev triggers
        # updating the plot
        self.ui.numPoints.returnPressed.connect(self.points_entered)
        self.ui.numStdDevs.returnPressed.connect(self.stdDevEntered)

        # Triggers a redrawing upon pressing enter in the search bar.
        # Proper usage should be using the search bar to search, and selecting
        # from the results in the list. If it's not in the list, it's an invalid
        # PV with no reason to attempt plotting
        self.ui.searchInputA.returnPressed.connect(self.inputActivated)
        self.ui.searchInputB.returnPressed.connect(self.inputActivated)

    def showGrid(self):
        self.plot.showGrid(self.ui.checkBoxShowGrid.isChecked(),
                           self.ui.checkBoxShowGrid.isChecked())

    def setUpGraph(self):
        layout = QGridLayout()
        self.ui.widgetPlot.setLayout(layout)
        layout.addWidget(self.plot, 0, 0)
        self.plot.showGrid(1, 1)

    def loadStyleSheet(self):
        currDir = path.abspath(path.dirname(__file__))
        cssFile = path.join(currDir, "style.css")
        try:
            with open(cssFile, "r") as f:
                self.setStyleSheet(f.read())
        except IOError:
            print("Error loading style sheet")
            pass

    def create_status_bar(self):
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
        self.search(self.ui.searchInputA, self.ui.bsaListA)

    def searchB(self):
        self.search(self.ui.searchInputB, self.ui.bsaListB)

    def setEnter(self, widget, enter, search, enter_rb):
        selection = widget.currentItem()
        enter.textChanged.disconnect()
        enter.setText(selection.text())
        QObject.connect(enter, SIGNAL("textChanged(const QString&)"), search)

        if not self.abort and enter_rb.isChecked():
            self.stop()
            self.timer.singleShot(250, self.initializePlot)

    def setEnterA(self):
        self.setEnter(self.ui.bsaListA, self.ui.searchInputA, self.searchA,
                      self.ui.searchButtonA)

    def setEnterB(self):
        self.setEnter(self.ui.bsaListB, self.ui.searchInputB, self.searchB,
                      self.ui.searchButtonB)

    def correctInput(self, errorMessage, acceptableTxt, textBox):
        self.statusBar().showMessage(errorMessage, 6000)
        textBox.setText(acceptableTxt)

    def correctNumpoints(self, errorMessage, acceptableValue):
        self.correctInput(errorMessage, str(acceptableValue), self.ui.numPoints)
        self.numPoints = acceptableValue

    def correctStdDevs(self, errorMessage, acceptableValue):
        self.correctInput(errorMessage, str(acceptableValue),
                          self.ui.numStdDevs)
        self.stdDevstoKeep = acceptableValue

    def stdDevEntered(self):
        try:
            self.stdDevstoKeep = float(self.ui.numStdDevs.text())
            if self.stdDevstoKeep <= 0:
                raise ValueError
        except ValueError:
            self.correctStdDevs('Enter a float > 0', 3.0)
            return

    def points_entered(self):
        try:
            self.numPoints = int(self.ui.numPoints.text())
        except ValueError:
            self.correctNumpoints('Enter an integer, 1 to 2800', 120)
            return

        if self.numPoints > 2800:
            self.correctNumpoints('Max # points is 2800', 2800)
            return

        if self.numPoints < 1:
            self.correctNumpoints('Min # points is 1', 1)
            return

        self.reinitialize_plot()

    ############################################################################
    # Where the magic happens (well, where it starts to happen). This
    # initializes the BSA plotting and then starts a timer to update the plot.
    ############################################################################
    def initializePlot(self):
        plotTypeIsValid = (self.ui.checkBoxAvsT.isChecked()
                           or self.ui.checkBoxBvsA.isChecked()
                           or self.ui.checkBoxFFT.isChecked())

        if not plotTypeIsValid:
            self.statusBar().showMessage('Pick a Plot Type (PV vs. time '
                                         'or B vs A)', 10000)
            return

        self.ui.startButton.setDisabled(True)
        self.abort = False

        self.cleanPlot()
        self.pvObjects["A"], self.pvObjects["B"] = None, None

        # Plot history buffer for one PV
        if self.ui.checkBoxAvsT.isChecked():
            if self.populateDevices(self.ui.dropdownButtonA, self.ui.dropdownA,
                                    self.ui.searchButtonA, self.ui.searchInputA, "A"):
                self.genPlotAndSetTimer(self.genTimePlotA,
                                        self.updateTimePlotA)

        # Plot for 2 PVs
        elif self.ui.checkBoxBvsA.isChecked():
            if self.updateValsFromInput():
                self.genPlotAndSetTimer(self.genPlotAB,
                                        self.updatePlotAB)

        # Plot power spectrum
        else:
            if self.populateDevices(self.ui.dropdownButtonA, self.ui.dropdownA,
                                    self.ui.searchButtonA, self.ui.searchInputA, "A"):
                self.genPlotAndSetTimer(self.InitializeFFTPlot,
                                        self.updatePlotFFT)

    def populateDevices(self, common_rb, common, enter_rb, enter, device):

        if common_rb.isChecked():
            self.devices[device] = str(common.currentText())

        elif enter_rb.isChecked():
            pv = str(enter.text()).strip()

            # Checks that it's non empty and that it's a BSA pv
            if pv and pv in self.bsapvs:
                self.devices[device] = pv
            else:
                self.printStatus('Device ' + device + ' invalid. Aborting.')
                self.ui.startButton.setEnabled(True)
                return False

        return True

    def printStatus(self, message, printToXterm=True):
        if printToXterm:
            print message

        self.statusBar().showMessage(message)

    def updateValsFromInput(self):

        if not self.populateDevices(self.ui.dropdownButtonA, self.ui.dropdownA,
                                    self.ui.searchButtonA, self.ui.searchInputA, "A"):
            return False

        if not self.populateDevices(self.ui.dropdownButtonB, self.ui.dropdownB,
                                    self.ui.searchButtonB, self.ui.searchInputB, "B"):
            return False

        self.printStatus("Initializing/Synchronizing " + self.devices["A"]
                         + " vs. " + self.devices["B"] + " buffers...")

        # Initial population of our buffers using the HSTBR PV's in our
        # callback functions
        self.clearAndUpdateCallbacks("HSTBR", resetTime=True)

        while ((not self.timeStamps["A"] or not self.timeStamps["B"])
               and not self.abort):
            QApplication.processEvents()

        self.adjustSynchronizedBuffers(True)

        # Switch to BR PVs to avoid pulling an entire history buffer on every
        # update.
        self.clearAndUpdateCallbacks("BR", resetRawBuffer=True)

        return True

    def clearAndUpdateCallbacks(self, suffix, resetTime=False,
                                resetRawBuffer=False):
        self.clearAndUpdateCallback("A", suffix, self.callbackA,
                                    self.devices["A"], resetTime,
                                    resetRawBuffer)
        self.clearAndUpdateCallback("B", suffix, self.callbackB,
                                    self.devices["B"], resetTime,
                                    resetRawBuffer)

    # noinspection PyTypeChecker
    def clearAndUpdateCallback(self, device, suffix, callback, pvName,
                               resetTime=False, resetRawBuffer=False):
        self.clearPV(device)

        # Without the time parameter, we wouldn't get the timestamp
        self.pvObjects[device] = PV(pvName + suffix, form='time')

        if resetTime:
            self.timeStamps[device] = None

        # Make sure that the initial raw buffer is synchronized and pad with
        # nans if it's less than 2800 points long
        if resetRawBuffer:
            nanArray = empty(2800 - self.synchronizedBuffers[device].size)
            nanArray[:] = nan
            self.rawBuffers[device] = \
                concatenate([self.synchronizedBuffers[device], nanArray])

        self.pvObjects[device].add_callback(callback)

    # Callback function for Device A
    # noinspection PyUnusedLocal
    def callbackA(self, pvname=None, value=None, timestamp=None, **kw):
        self.updateTimeAndBuffer("A", pvname, timestamp, value)

    # Callback function for Device B
    # noinspection PyUnusedLocal
    def callbackB(self, pvname=None, value=None, timestamp=None, **kw):
        self.updateTimeAndBuffer("B", pvname, timestamp, value)

    ############################################################################
    # This is where the data is actually acquired and saved to the buffers.
    # Callbacks are effectively listeners that listen for change, so we
    # basically put a callback on the PVs of interest (devices A and/or B) so
    # that every time the value of that PV changes, we get that new value and
    # append it to our raw data buffer for that device.
    # Initialization of the buffer is slightly different in that the listener is
    # put on the history buffer of that PV (denoted by the HSTBR suffix), so
    # that we just immediately write the previous 2800 points to our raw buffer
    ############################################################################
    def updateTimeAndBuffer(self, device, pvname, timestamp, value):

        def padRawBufferWithNans(start, end):
            rtbsaUtils.padWithNans(self.rawBuffers[device], start, end)

        if "HSTBR" in pvname:
            self.timeStamps[device] = timestamp

            # value is the buffer because we're monitoring the HSTBR PV
            self.rawBuffers[device] = value

            # Reset the counter every time we reinitialize the plot
            self.counter[device] = 0

        else:
            if not self.timeStamps[device]:
                return

            rate = self.getRate()
            if rate < 1:
                return

            scalingFactor = 1.0 / rate
            elapsedPulses = round((timestamp - self.timeStamps[device])
                                  / scalingFactor)

            currIdx = int((timestamp / scalingFactor) % self.numPoints)

            if elapsedPulses <= 0:
                return

            # Pad the buffer with nans for missed pulses
            elif elapsedPulses > 1:

                # noinspection PyTypeChecker
                lastIdx = int((self.timeStamps[device] / scalingFactor)
                              % self.numPoints)

                # Take care of wrap around
                if currIdx < lastIdx:
                    padRawBufferWithNans(lastIdx + 1, self.numPoints)
                    padRawBufferWithNans(0, currIdx)

                else:
                    padRawBufferWithNans(lastIdx + 1, currIdx)

            # Directly index into the raw buffer using the timestamp
            self.rawBuffers[device][currIdx] = value

            self.counter[device] += elapsedPulses
            self.timeStamps[device] = timestamp
            self.currIdx[device] = currIdx

    def clearPV(self, device):
        pv = self.pvObjects[device]
        if pv:
            pv.clear_callbacks()
            pv.disconnect()

    def adjustSynchronizedBuffers(self, syncByTime=False):
        numBadShots = self.populateSynchronizedBuffers(syncByTime)
        blength = 2800 - numBadShots

        # Make sure the buffer size doesn't exceed the desired number of points
        if self.numPoints < blength:
            self.synchronizedBuffers["A"] = \
                self.synchronizedBuffers["A"][:self.numPoints]

            self.synchronizedBuffers["B"] = \
                self.synchronizedBuffers["B"][:self.numPoints]

    # A spin loop that waits until the beam rate is at least 1Hz
    def waitForRate(self):

        start_time = time()
        gotStuckAndNeedToUpdateMessage = False

        # self.rate is a PV, such that .value is shorthand for .getval
        while self.ratePV.value < 2:
            # noinspection PyArgumentList
            QApplication.processEvents()

            if time() - start_time > 0.5:
                gotStuckAndNeedToUpdateMessage = True
                self.printStatus("Waiting for beam rate to be at least 1Hz...",
                                 False)

        if gotStuckAndNeedToUpdateMessage:
            self.printStatus("Running", False)

        return rtbsaUtils.rateDict[self.ratePV.value]

    ############################################################################
    # Time 1 is when Device A started acquiring data, and Time 2 is when Device
    # B started acquiring data. Since they're not guaranteed to start
    # acquisition at the same time, one data buffer might be ahead of the other,
    # meaning that the intersection of the two buffers would not include the
    # first n elements of one and the last n elements of the other. See the
    # diagram below, where the dotted line represents the time axis (one buffer
    # is contained  by square brackets [], the other by curly braces {}, and the
    # times where each starts and ends is indicated right underneath).
    #
    #
    #          [           {                            ]           }
    # <----------------------------------------------------------------------> t
    #       t1_start    t2_start                     t1_end      t2_end
    #
    #
    # Note that both buffers are of the same size (2800) so that:
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
    #
    # That whole rigmarole only applies to the initial population of the buffers
    # (where we're pulling the entire history buffer at once using the HSTBR
    # suffix). From then on, we're indexing into the raw buffers using the
    # pulse ID modulo 2800, so they're inherently synchronized
    ############################################################################
    def populateSynchronizedBuffers(self, syncByTime):

        def padSyncBufferWithNans(device, startIdx, endIdx):
            rtbsaUtils.padWithNans(self.synchronizedBuffers[device],
                                   startIdx - 2, endIdx + 1)

        def checkIndices(device, startIdx, endIdx):
            # Check for index wraparound
            if endIdx < startIdx:
                padSyncBufferWithNans(device, startIdx, self.numPoints - 1)
                padSyncBufferWithNans(device, 0, endIdx)
            else:
                padSyncBufferWithNans(device, startIdx, endIdx)

        if syncByTime:
            numBadShots = int(round((self.timeStamps["B"] - self.timeStamps["A"])
                                    * self.getRate()))

            startA, endA = rtbsaUtils.getIndices(numBadShots, 1)
            startB, endB = rtbsaUtils.getIndices(numBadShots, -1)

            self.synchronizedBuffers["A"] = self.rawBuffers["A"][startA:endA]
            self.synchronizedBuffers["B"] = self.rawBuffers["B"][startB:endB]

            return abs(numBadShots)
        else:

            self.synchronizedBuffers["A"] = self.rawBuffers["A"]
            self.synchronizedBuffers["B"] = self.rawBuffers["B"]

            # The timestamps and indices get updated by the callbacks, so we
            # store the values at the time of buffer-copying
            timeStampA = self.timeStamps["A"]
            timeStampB = self.timeStamps["B"]

            currIdxA = self.currIdx["A"]
            currIdxB = self.currIdx["B"]

            diff = timeStampA - timeStampB

            if diff > 0:
                checkIndices("A", currIdxB, currIdxA)

            elif diff < 0:
                checkIndices("B", currIdxA, currIdxB)

            return 0

    def genPlotAndSetTimer(self, genPlot, updateMethod):
        if self.abort:
            return

        try:
            genPlot()
        except UnboundLocalError:
            self.printStatus('No Data, Aborting Plotting Algorithm')
            return

        self.timer = QTimer(self)

        # Run updateMethod every updatetime milliseconds
        self.timer.singleShot(self.updateTime, updateMethod)

        self.printStatus('Running')

    # noinspection PyTypeChecker
    def genTimePlotA(self):
        newData = self.initializeData()

        if not newData.size:
            self.printStatus('Invalid PV? Unable to get data. Aborting.')
            self.ui.startButton.setEnabled(True)
            return

        data = newData[:self.numPoints]

        self.plotAttributes["curve"] = PlotCurveItem(data, pen=1)
        self.plot.addItem(self.plotAttributes["curve"])

        self.plotFit(arange(self.numPoints), data, self.devices["A"])

    ############################################################################
    # This is the main plotting function for "Plot A vs Time" that gets called
    # every self.updateTime seconds
    # noinspection PyTypeChecker
    ############################################################################
    def updateTimePlotA(self):

        if not self.checkPlotStatus():
            return

        xData, yData = self.filterTimePlotBuffer()

        if yData.size:
            self.plotAttributes["curve"].setData(yData)
            if self.ui.checkBoxAutoscale.isChecked():
                mx = max(yData)
                mn = min(yData)
                if mx - mn > .00001:
                    self.plot.setYRange(mn, mx)
                    self.plot.setXRange(0, len(yData))

            if self.ui.checkBoxShowAve.isChecked():
                rtbsaUtils.setPosAndText(self.text["avg"], mean(yData), 0,
                                         min(yData), 'AVG: ')

            if self.ui.checkBoxShowStdDev.isChecked():
                rtbsaUtils.setPosAndText(self.text["std"], std(yData),
                                         self.numPoints / 4, min(yData),
                                         'STD: ')

            if self.ui.checkBoxCorrCoeff.isChecked():
                self.text["corr"].setText('')

            if self.ui.checkBoxLinFit.isChecked():
                self.text["slope"].setPos(self.numPoints / 2, min(yData))
                self.getLinearFit(xData, yData, True)

            elif self.ui.checkBoxPolyFit.isChecked():
                self.text["slope"].setPos(self.numPoints / 2, min(yData))
                self.getPolynomialFit(xData, yData, True)

        self.timer.singleShot(self.updateTime, self.updateTimePlotA)

    def checkPlotStatus(self):
        QApplication.processEvents()

        if self.abort:
            return False

        self.waitForRate()

        # kill switch to stop backgrounded, forgetten GUIs. Somewhere in the
        # ballpark of 20 minutes assuming 120Hz
        if self.counter["A"] > 150000:
            self.stop()
            self.printStatus("Stopping due to inactivity")

        return True

    def filterTimePlotBuffer(self):
        currIdx = self.currIdx["A"]
        choppedBuffer = self.rawBuffers["A"][:self.numPoints]

        # All this nonsense to make it scroll :P Thanks to Ben for the
        # inspiration!
        if currIdx > 0:
            choppedBuffer = concatenate([choppedBuffer[currIdx:],
                                         choppedBuffer[:currIdx]])

        xData, yData = rtbsaUtils.filterBuffers(choppedBuffer,
                                                lambda x: ~isnan(x),
                                                arange(self.numPoints),
                                                choppedBuffer)

        if self.devices["A"] == "BLEN:LI24:886:BIMAX":
            xData, yData = rtbsaUtils.filterBuffers(yData,
                                                    lambda x:
                                                    x < rtbsaUtils.IPK_LIMIT,
                                                    xData, yData)

        if self.ui.checkBoxStdDev.isChecked():
            stdDevFilterFunc = self.StdDevFilterFunc(mean(yData), std(yData))
            xData, yData = rtbsaUtils.filterBuffers(yData, stdDevFilterFunc,
                                                    xData, yData)
        return xData, yData

    def getLinearFit(self, xData, yData, updateExistingPlot):
        # noinspection PyTupleAssignmentBalance
        m, b = polyfit(xData, yData, 1)
        fitData = polyval([m, b], xData)

        self.text["slope"].setText('Slope: ' + str("{:.3e}".format(m)))

        if updateExistingPlot:
            self.plotAttributes["fit"].setData(xData, fitData)
        else:
            # noinspection PyTypeChecker
            self.plotAttributes["fit"] = PlotCurveItem(xData, fitData, 'g-',
                                                       linewidth=1)

    def getPolynomialFit(self, xData, yData, updateExistingPlot):
        co = polyfit(xData, yData, self.fitOrder)
        pol = poly1d(co)
        xDataSorted = sorted(xData)
        fit = pol(xDataSorted)

        if updateExistingPlot:
            self.plotAttributes["parab"].setData(xDataSorted, fit)
        else:
            # noinspection PyTypeChecker
            self.plotAttributes["parab"] = PlotCurveItem(xDataSorted, fit,
                                                         pen=3, size=2)

        if self.fitOrder == 2:
            self.text["slope"].setText('Peak: ' + str(-co[1] / (2 * co[0])))

        elif self.fitOrder == 3:
            self.text["slope"].setText(str("{:.2e}".format(co[0])) + 'x^3'
                                       + str("+{:.2e}".format(co[1]))
                                       + 'x^2'
                                       + str("+{:.2e}".format(co[2])) + 'x'
                                       + str("+{:.2e}".format(co[3])))


    def genPlotAB(self):
        if self.ui.checkBoxStdDev.isChecked():
            self.plotCurveAndFit(self.filteredBuffers["A"],
                                 self.filteredBuffers["B"])
        else:
            self.plotCurveAndFit(self.synchronizedBuffers["A"],
                                 self.synchronizedBuffers["B"])

    def plotCurveAndFit(self, xData, yData):
        # noinspection PyTypeChecker
        self.plotAttributes["curve"] = ScatterPlotItem(xData, yData, pen=1,
                                                       symbol='x', size=5)
        self.plot.addItem(self.plotAttributes["curve"])
        self.plotFit(xData, yData,
                     self.devices["B"] + ' vs. ' + self.devices["A"])

    def plotFit(self, xData, yData, title):
        self.plot.addItem(self.plotAttributes["curve"])
        self.plot.setTitle(title)

        # Fit line
        if self.ui.checkBoxLinFit.isChecked():
            self.getLinearFit(xData, yData, False)
            self.plot.addItem(self.plotAttributes["fit"])

        # Fit polynomial
        elif self.ui.checkBoxPolyFit.isChecked():
            self.ui.fitOrder.setDisabled(False)
            try:
                self.getPolynomialFit(xData, yData, False)
                self.plot.addItem(self.plotAttributes["parab"])
            except linalg.linalg.LinAlgError:
                print "Error getting polynomial fit"

    ############################################################################
    # This is the main plotting function for "Plot B vs A" that gets called
    # every self.updateTime milliseconds
    ############################################################################
    def updatePlotAB(self):
        if not self.checkPlotStatus():
            return

        QApplication.processEvents()

        self.adjustSynchronizedBuffers()
        self.filterNans()
        self.filterPeakCurrent()

        if self.ui.checkBoxStdDev.isChecked():
            self.filterStdDev()
            self.updateLabelsAndFit(self.filteredBuffers["A"],
                                    self.filteredBuffers["B"])
        else:
            self.updateLabelsAndFit(self.synchronizedBuffers["A"],
                                    self.synchronizedBuffers["B"])

        self.timer.singleShot(self.updateTime, self.updatePlotAB)

    def filterNans(self):
        def filterFunc(x): return ~isnan(x)

        self.filterData(self.synchronizedBuffers["A"], filterFunc, True)
        self.filterData(self.synchronizedBuffers["B"], filterFunc, True)

    # Need to filter out errant indices from both buffers to keep them
    # synchronized
    def filterData(self, dataBuffer, filterFunc, changeOriginal):
        bufferA, bufferB = rtbsaUtils.filterBuffers(dataBuffer, filterFunc,
                                                    self.synchronizedBuffers[
                                                        "A"],
                                                    self.synchronizedBuffers[
                                                        "B"])

        if changeOriginal:
            self.synchronizedBuffers["A"] = bufferA
            self.synchronizedBuffers["B"] = bufferB
        else:
            self.filteredBuffers["A"] = bufferA
            self.filteredBuffers["B"] = bufferB

    # This PV gets insane values, apparently
    def filterPeakCurrent(self):
        def filterFunc(x):
            return x < rtbsaUtils.IPK_LIMIT

        if self.devices["A"] == "BLEN:LI24:886:BIMAX":
            self.filterData(self.synchronizedBuffers["A"], filterFunc, True)
        if self.devices["B"] == "BLEN:LI24:886:BIMAX":
            self.filterData(self.synchronizedBuffers["B"], filterFunc, True)

    def filterStdDev(self):

        bufferA = self.synchronizedBuffers["A"]
        bufferB = self.synchronizedBuffers["B"]

        self.filterData(bufferA, self.StdDevFilterFunc(mean(bufferA),
                                                       std(bufferA)), False)

        self.filterData(bufferB, self.StdDevFilterFunc(mean(bufferB),
                                                       std(bufferB)), False)

    def StdDevFilterFunc(self, average, stdDev):
        return lambda x: abs(x - average) < self.stdDevstoKeep * stdDev

    # noinspection PyTypeChecker
    def updateLabelsAndFit(self, bufferA, bufferB):
        self.plotAttributes["curve"].setData(bufferA, bufferB)

        try:
            if self.ui.checkBoxAutoscale.isChecked():
                self.setPlotRanges(bufferA, bufferB)

            minBufferA = nanmin(bufferA)
            minBufferB = nanmin(bufferB)
            maxBufferA = nanmax(bufferA)
            maxBufferB = nanmax(bufferB)

            if self.ui.checkBoxShowAve.isChecked():
                rtbsaUtils.setPosAndText(self.text["avg"], nanmean(bufferB),
                                         minBufferA,
                                         minBufferB, 'AVG: ')

            if self.ui.checkBoxShowStdDev.isChecked():
                xPos = (minBufferA + (minBufferA + maxBufferA) / 2) / 2

                rtbsaUtils.setPosAndText(self.text["std"], nanstd(bufferB),
                                         xPos, minBufferB, 'STD: ')

            if self.ui.checkBoxCorrCoeff.isChecked():
                correlation = corrcoef(bufferA, bufferB)
                rtbsaUtils.setPosAndText(self.text["corr"], correlation.item(1),
                                         minBufferA, maxBufferB,
                                         "Corr. Coefficient: ")

            if self.ui.checkBoxLinFit.isChecked():
                self.text["slope"].setPos((minBufferA + maxBufferA) / 2,
                                          minBufferB)
                self.getLinearFit(bufferA, bufferB, True)

            elif self.ui.checkBoxPolyFit.isChecked():
                self.text["slope"].setPos((minBufferA + maxBufferA) / 2,
                                          minBufferB)
                self.getPolynomialFit(bufferA, bufferB, True)

        except ValueError:
            print "Error updating plot range"

    def setPlotRanges(self, bufferA, bufferB):
        mx = nanmax(bufferB)
        mn = nanmin(bufferB)

        if mn != mx:
            self.plot.setYRange(mn, mx)

        mx = nanmax(bufferA)
        mn = nanmin(bufferA)

        if mn != mx:
            self.plot.setXRange(mn, mx)

    def InitializeFFTPlot(self):
        self.genPlotFFT(self.initializeData(), False)

    # TODO I have no idea what's happening here
    def genPlotFFT(self, newdata, updateExistingPlot):

        if not newdata.size:
            return None

        newdata = newdata[:self.numPoints]

        nans, x = isnan(newdata), lambda z: z.nonzero()[0]
        # interpolate nans
        newdata[nans] = interp(x(nans), x(~nans), newdata[~nans])
        # remove DC component
        newdata = newdata - mean(newdata)

        newdata = concatenate([newdata, zeros(self.numPoints * 2)])

        ps = abs(fft.fft(newdata)) / newdata.size

        self.waitForRate()
        frequencies = fft.fftfreq(newdata.size, 1.0 / self.getRate())
        keep = (frequencies >= 0)
        ps = ps[keep]
        frequencies = frequencies[keep]
        idx = argsort(frequencies)

        if updateExistingPlot:
            self.plotAttributes["curve"].setData(x=frequencies[idx], y=ps[idx])
        else:
            # noinspection PyTypeChecker
            self.plotAttributes["curve"] = PlotCurveItem(x=frequencies[idx],
                                                         y=ps[idx], pen=1)

        self.plot.addItem(self.plotAttributes["curve"])
        self.plot.setTitle(self.devices["A"])
        self.plotAttributes["frequencies"] = frequencies

        return ps

    # noinspection PyTypeChecker
    def cleanPlot(self):
        self.plot.clear()

        self.text["avg"] = TextItem('', color=(200, 200, 250), anchor=(0, 1))
        self.text["std"] = TextItem('', color=(200, 200, 250), anchor=(0, 1))
        self.text["slope"] = TextItem('', color=(200, 200, 250), anchor=(0, 1))
        self.text["corr"] = TextItem('', color=(200, 200, 250), anchor=(0, 1))

        plotLabels = [self.text["avg"], self.text["std"], self.text["slope"],
                      self.text["corr"]]

        for plotLabel in plotLabels:
            self.plot.addItem(plotLabel)

    def initializeData(self):
        self.printStatus("Initializing " + self.devices["A"] + " buffer...",
                         True)

        if self.ui.dropdownButtonA.isChecked():
            self.devices["A"] = str(self.ui.dropdownA.currentText())

        elif self.ui.searchButtonA.isChecked():
            pv = str(self.ui.searchInputA.text()).strip()
            if pv and pv in self.bsapvs:
                self.devices["A"] = pv
            else:
                return None
        else:
            return None

        # Initializing our data by putting a callback on the history buffer PV
        self.clearAndUpdateCallback("A", "HSTBR", self.callbackA,
                                    self.devices["A"], True)

        while (not self.timeStamps["A"]) and not self.abort:
            QApplication.processEvents()

        # Removing that callback and manually appending new values to our local
        # data buffer using the usual PV
        # TODO ask Ahmed what the BR is for
        self.clearAndUpdateCallback("A", "BR", self.callbackA,
                                    self.devices["A"])

        # This was populated in the callback function
        return self.rawBuffers["A"]

    ############################################################################
    # This is the main plotting function for "Plot A FFT" that gets called
    # every self.updateTime seconds
    ############################################################################
    def updatePlotFFT(self):
        if not self.checkPlotStatus():
            return

        ps = self.genPlotFFT(self.rawBuffers["A"], True)

        if self.ui.checkBoxAutoscale.isChecked():
            mx = max(ps)
            mn = min(ps)
            if mx - mn > .00001:
                frequencies = self.plotAttributes["frequencies"]
                self.plot.setYRange(mn, mx)
                # noinspection PyTypeChecker
                self.plot.setXRange(min(frequencies), max(frequencies))

        self.timer.singleShot(self.updateTime, self.updatePlotFFT)

    def AvsTClick(self):
        if not self.ui.checkBoxAvsT.isChecked():
            pass
        else:
            self.ui.checkBoxBvsA.setChecked(False)
            self.ui.checkBoxFFT.setChecked(False)
            self.AvsBClick()

    def AvsBClick(self):
        if not self.ui.checkBoxBvsA.isChecked():
            self.ui.groupBoxB.setDisabled(True)
            self.ui.searchButtonB.setChecked(False)
            self.ui.searchButtonB.setDisabled(True)
            self.ui.searchInputB.setDisabled(True)
            self.ui.dropdownB.setDisabled(True)
            self.ui.dropdownButtonB.setChecked(False)
            self.ui.dropdownButtonB.setDisabled(True)
        else:
            self.ui.checkBoxAvsT.setChecked(False)
            self.ui.checkBoxFFT.setChecked(False)
            self.AvsTClick()
            self.ui.groupBoxB.setDisabled(False)
            self.ui.bsaListB.setDisabled(True)
            self.ui.searchButtonB.setDisabled(False)
            self.ui.searchInputB.setDisabled(True)
            self.ui.dropdownButtonB.setDisabled(False)
            self.ui.dropdownButtonB.setChecked(True)
            self.ui.dropdownB.setDisabled(False)

        self.stop()
        self.timer.singleShot(250, self.initializePlot)

    def AFFTClick(self):
        if not self.ui.checkBoxFFT.isChecked():
            pass
        else:
            self.ui.checkBoxBvsA.setChecked(False)
            self.ui.checkBoxAvsT.setChecked(False)
            self.AvsBClick()

    def avg_click(self):
        if not self.ui.checkBoxShowAve.isChecked():
            self.text["avg"].setText('')

    def std_click(self):
        if not self.ui.checkBoxShowStdDev.isChecked():
            self.text["std"].setText('')

    def corr_click(self):
        if not self.ui.checkBoxCorrCoeff.isChecked():
            self.text["corr"].setText('')

    def enter_1_click(self):
        if self.ui.searchButtonA.isChecked():
            self.ui.searchInputA.setDisabled(False)
            self.ui.bsaListA.setDisabled(False)
            self.ui.dropdownButtonA.setChecked(False)
            self.ui.dropdownA.setDisabled(True)
        else:
            self.ui.searchInputA.setDisabled(True)

    def enter_2_click(self):
        if self.ui.searchButtonB.isChecked():
            self.ui.searchInputB.setDisabled(False)
            self.ui.bsaListB.setDisabled(False)
            self.ui.dropdownButtonB.setChecked(False)
            self.ui.dropdownB.setDisabled(True)
        else:
            self.ui.searchInputB.setDisabled(True)

    def common_1_click(self):
        if self.ui.dropdownButtonA.isChecked():
            self.ui.dropdownA.setEnabled(True)
            self.ui.searchButtonA.setChecked(False)
            self.ui.searchInputA.setDisabled(True)
            self.ui.bsaListA.setDisabled(True)
        else:
            self.ui.dropdownA.setEnabled(False)
        self.inputActivated()

    def inputActivated(self):
        if not self.abort:
            self.stop()
            self.timer.singleShot(250, self.initializePlot)

    def common_2_click(self):
        if self.ui.dropdownButtonB.isChecked():
            self.ui.dropdownB.setEnabled(True)
            self.ui.searchButtonB.setChecked(False)
            self.ui.searchInputB.setDisabled(True)
            self.ui.bsaListB.setDisabled(True)
        else:
            self.ui.dropdownB.setEnabled(False)
        self.inputActivated()

    def line_click(self):
        self.ui.checkBoxPolyFit.setChecked(False)
        self.ui.fitOrder.setDisabled(True)
        self.ui.labelFitOrder.setDisabled(True)
        self.reinitialize_plot()

    def fitOrderActivated(self):
        try:
            self.fitOrder = int(self.ui.fitOrder.text())
        except ValueError:
            self.statusBar().showMessage('Enter an integer, 1-10', 6000)
            return

        if self.fitOrder > 10 or self.fitOrder < 1:
            self.statusBar().showMessage('Really?  That is going to be useful'
                                         + ' to you?  The (already ridiculous)'
                                         + ' range is 1-10.  Hope you win a '
                                         + 'nobel prize jackass.', 6000)
            self.ui.fitOrder.setText('2')
            self.fitOrder = 2

        if self.fitOrder != 2:
            try:
                self.text["slope"].setText('')
            except AttributeError:
                pass

    def parab_click(self):
        self.ui.checkBoxLinFit.setChecked(False)

        if not self.ui.checkBoxPolyFit.isChecked():
            self.ui.fitOrder.setDisabled(True)
            self.ui.labelFitOrder.setDisabled(True)
        else:
            self.ui.fitOrder.setEnabled(True)
            self.ui.labelFitOrder.setEnabled(True)

        self.reinitialize_plot()

    # This is a mess, but it works (used if user changes number points,
    # fit type etc.)
    def reinitialize_plot(self):
        self.cleanPlot()

        # Setup for single PV plotting
        if self.ui.checkBoxAvsT.isChecked():
            self.genTimePlotA()

        elif self.ui.checkBoxBvsA.isChecked():
            self.genPlotAB()
        else:
            self.genPlotFFT(self.synchronizedBuffers["A"], False)

    def logbook(self):
        rtbsaUtils.logbook('Python Real-Time BSA', 'BSA Data',
                           str(self.numPoints) + ' points', self.plot.plotItem)
        self.statusBar().showMessage('Sent to LCLS Physics Logbook!', 10000)

    def MCCLog(self):
        rtbsaUtils.MCCLog('/tmp/RTBSA.png', '/tmp/RTBSA.ps', self.plot.plotItem)

    def clearCallbacks(self, device):
        if self.pvObjects[device]:
            self.pvObjects[device].clear_callbacks()
            self.pvObjects[device].disconnect()

    def stop(self):
        self.clearCallbacks("A")

        if self.pvObjects["B"]:
            self.clearCallbacks("B")

        self.abort = True
        self.statusBar().showMessage('Stopped')
        self.ui.startButton.setDisabled(False)
        QApplication.processEvents()

    def create_menu(self):

        load_file_action = self.create_action("&Save plot", shortcut="Ctrl+S",
                                              slot=self.save_plot,
                                              tip="Save the plot")

        quit_action = self.create_action("&Quit", slot=self.close,
                                         shortcut="Ctrl+Q",
                                         tip="Close the application")

        rtbsaUtils.add_actions(self.file_menu, (load_file_action, None,
                                                quit_action))

        about_action = self.create_action("&About", shortcut='F1',
                                          slot=self.on_about, tip='About')

        rtbsaUtils.add_actions(self.help_menu, (about_action,))

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
        # noinspection PyTypeChecker,PyCallByClass
        filePath = unicode(QFileDialog.getSaveFileName(self, 'Save file', '',
                                                       file_choices))
        if filePath:
            self.ui.widgetPlot.canvas.print_figure(filePath, dpi=100)
            self.statusBar().showMessage('Saved to %s' % filePath, 2000)

    def on_about(self):
        msg = ("Can you read this?  If so, congratulations. You are a magical, "
               + "marvelous troll.")
        # noinspection PyCallByClass
        QMessageBox.about(self, "About", msg.strip())


# TODO I bless the rains down in Africa!
def main():
    app = QApplication(argv)
    window = RTBSA()
    window.show()
    exit(app.exec_())


if __name__ == "__main__":
    main()
