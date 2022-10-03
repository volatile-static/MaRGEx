"""
Created on Thu June 2 2022
@author: J.M. Algarín, MRILab, i3M, CSIC, Valencia
@email: josalggui@i3m.upv.es
@Summary: code to obtain a good combination of tuning/matching
Specific hardware from MRILab @ i3M is required
"""

import os
import sys
#*****************************************************************************
# Add path to the working directory
path = os.path.realpath(__file__)
ii = 0
for char in path:
    if (char=='\\' or char=='/') and path[ii+1:ii+14]=='PhysioMRI_GUI':
        sys.path.append(path[0:ii+1]+'PhysioMRI_GUI')
        sys.path.append(path[0:ii+1]+'marcos_client')
    ii += 1
#******************************************************************************
from PyQt5.QtCore import QThreadPool
import experiment as ex
import numpy as np
import seq.mriBlankSeq as blankSeq  # Import the mriBlankSequence for any new sequence.
from worker import Worker
import serial
import time
import configs.hw_config as hw

from plotview.spectrumplot import SpectrumPlot
import pyqtgraph as pg


class AutoTuning(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(AutoTuning, self).__init__()
        # Input the parameters
        self.expt = None
        self.repeat = None
        self.threadpool = None
        self.rfExAmp = None
        self.rfExTime = None
        self.txChannel = None
        self.freqOffset = None
        self.addParameter(key='seqName', string='AutoTuningInfo', val='AutoTuning')
        self.addParameter(key='freqOffset', string='RF frequency offset (kHz)', val=0.0, field='RF')
        self.addParameter(key='rfExTime', string='RF excitation time (s)', val=10.0, field='RF')
        self.addParameter(key='rfExAmp', string='RF excitation amplitude (a.u.)', val=0.02, field='RF')
        self.addParameter(key='txChannel', string='Tx channel', val=0, field='RF')

    def sequenceInfo(self):
        print("\n RF Auto-tuning")
        print("Author: Dr. J.M. Algarín")
        print("Contact: josalggui@i3m.upv.es")
        print("mriLab @ i3M, CSIC, Spain")
        print("Look for the best combination of tuning/matching.")
        print("Specific hardware from MRILab @ i3M is required. \n")

    def sequenceTime(self):
        return 0  # minutes, scanTime

    def sequenceRun(self, plotSeq=0):
        # Create the inputs automatically as class properties
        for key in self.mapKeys:
            setattr(self, key, self.mapVals[key])

        # Fix units to MHz and us
        hw.larmorFreq = 3.066 # MHz
        self.freqOffset *= 1e-3  # MHz
        self.rfExTime *= 1e6 # us

        # # SEQUENCE
        self.expt = ex.Experiment(lo_freq=hw.larmorFreq + self.freqOffset, init_gpa=False)
        t0 = 5
        self.iniSequence(t0, np.array([0, 0, 0]))
        t0 = 10
        # self.rfRecPulse(t0, self.rfExTime, self.rfExAmp, txChannel=0)
        self.ttl(t0, self.rfExTime + hw.blkTime, channel=1)
        self.rfRawPulse(t0 + hw.blkTime, self.rfExTime, self.rfExAmp, txChannel=1)
        t0 += hw.blkTime + self.rfExTime + 10000
        self.endSequence(t0)

        # Run sequence continuously
        if not plotSeq:
            self.repeat = True
            # Sweep the tuning matching states in parallel thread
            self.threadpool = QThreadPool()
            print("Multithreading with maximum %d threads \n" % self.threadpool.maxThreadCount())
            worker = Worker(self.runAutoTuning)  # Any other args, kwargs are passed to the run function
            self.threadpool.start(worker)
            # Excite
            n = 0
            while self.repeat:
                print("Running...")
                self.expt.run()
                n += 1
                time.sleep(1)

            print("Ready!")
        self.expt.__del__()

    def sequenceAnalysis(self, obj=''):
        self.mapVals['bestSState'] = self.bestSState
        self.mapVals['bestTmState'] = self.bestTmState
        self.mapVals['minVoltage'] = self.vMin
        self.saveRawData()

        return([])

    def runTest(self):
        time.sleep(10)
        print("Soy A")
        self.repeat = False

    def runAutoTuning(self):
        start = time.time()

        # Open arduino serial port
        arduino = serial.Serial(port='COM7', baudrate=115200, timeout=.1)
        print('\n Arduino connected!')
        time.sleep(1)

        # Start states sweep
        arduino.write(b'.')

        # Read data
        while arduino.in_waiting == 0:
            time.sleep(0.1)
        time.sleep(0.2)
        result = arduino.readline().decode('utf-8')[0:-2].split(",")

        # Get the best state
        self.bestSState = result[0]
        self.bestTmState = result[1]
        self.vMin = int(result[2])/1023*5

        # Complete serial and tuning/matching binary values
        while len(self.bestSState)<5:
            self.bestSState += "0"
        while len(self.bestTmState)<10:
            self.bestTmState = "0"+self.bestTmState

        # Print the best state
        print("Best series state = ", self.bestSState)
        print("Best tuning/matching state = ", self.bestTmState)
        print("Min voltage = %0.3f V" % self.vMin)

        # Close and destroy the arduino
        arduino.close()
        arduino.__del__()

        # Switch the repeat variable
        self.repeat = False

        stop = time.time()
        print("Elapsed time = %0.1f s" % (stop-start))


if __name__ == '__main__':
    seq = AutoTuning()
    seq.sequenceRun()
    seq.sequenceAnalysis(obj='Standalone')
