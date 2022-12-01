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
from nanoVNA.Hardware import get_interfaces, get_VNA
import numpy as np
import seq.mriBlankSeq as blankSeq  # Import the mriBlankSequence for any new sequence.
from worker import Worker
import serial
import time
import configs.hw_config as hw

class AutoTuning(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(AutoTuning, self).__init__()
        # Input the parameters
        self.frequencies = None
        self.expt = None
        self.repeat = None
        self.threadpool = None
        self.rfExAmp = None
        self.rfExTime = None
        self.txChannel = None
        self.freqOffset = None
        self.addParameter(key='seqName', string='AutoTuningInfo', val='AutoTuning')
        self.addParameter(key='accuracy', string='Accuracy (dB)', val=-20.0, field='RF')
        self.addParameter(key='iterations', string='Max iterations', val=10, field='RF')

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

        # Run sequence continuously
        self.threadpool = QThreadPool()
        print("Multithreading with maximum %d threads \n" % self.threadpool.maxThreadCount())
        worker = Worker(self.runAutoTuning)  # Any other args, kwargs are passed to the run function
        self.threadpool.start(worker)

    def sequenceAnalysis(self, obj=''):
        # self.mapVals['bestSState'] = self.bestSState
        # self.mapVals['bestTmState'] = self.bestTmState
        # self.mapVals['minVoltage'] = self.vMin
        self.saveRawData()

        return([])

    def runTest(self):
        time.sleep(10)
        print("Soy A")
        self.repeat = False

    def runAutoTuning(self):
        start = time.time()
        nCap = 5

        # Combinations
        states = [None]*2**nCap
        zero = "00000"
        for state in range(2**nCap):
            prov = bin(state)[2::]
            x = len(prov)
            if x < nCap:
                states[state] = zero[0:nCap-x]+prov
            else:
                states[state] = prov

        # Series reactance
        cs = np.array([np.Inf, 8, 3.9, 1.8, 1])*1e-9
        statesCs = np.zeros(2**nCap)
        for state in range(2**nCap):
            for c in range(nCap):
                if c == 0 and int(states[state][0]) == 0:
                    statesCs[state] += 0
                else:
                    statesCs[state] += int(states[state][c])*cs[c]
        statesXs = -1/(2*np.pi*hw.larmorFreq*1e6*statesCs)

        # Tuning capacitor
        ct = np.array([326, 174, 87, 44, 26])*1e-12
        statesCt = np.zeros(2**nCap)
        for state in range(2**nCap):
            for c in range(nCap):
                statesCt[state] += int(states[state][c])*ct[c]

        # Matching capacitors
        cm = np.array([np.Inf, 500, 262, 142, 75])*1e-12
        statesCm = np.zeros(2 ** nCap)
        for state in range(2 ** nCap):
            for c in range(nCap):
                if c == 0 and int(states[state][0]) == 0:
                    statesCm[state] += 0
                else:
                    statesCm[state] += int(states[state][c]) * cm[c]
        statesXm = -1 / (2 * np.pi * hw.larmorFreq * 1e6 * statesCm)

        # Open nanoVNA
        # scan serial ports and connect
        interface = get_interfaces()[0]
        interface.open()
        interface.timeout = 0.05
        time.sleep(0.1)
        vna = get_VNA(interface)

        # # Open arduino serial port
        # arduino = serial.Serial(port='COM7', baudrate=115200, timeout=.1)
        # print('\n Arduino connected!')
        # time.sleep(1)
        #
        # # Initial state
        # arduino.write(b'100000000010000')
        #
        # # Read arduino state
        # while arduino.in_waiting == 0:
        #     time.sleep(0.1)
        # result = arduino.readline()

        # Measure the initial impedance
        if self.frequencies==None:
            self.frequencies = np.array(vna.readFrequencies())*1e-6
        idf = np.argmin(abs(self.frequencies-hw.larmorFreq))
        s11 = np.array([float(value) for value in vna.readValues("data 0")[idf].split(" ")])  # "data 0"->S11, "data 1"->S21
        s11 = s11[0]+s11[1]*1j
        impedance = 50*(1.+s11)/(1.-s11)
        r0 = impedance.real
        x0 = impedance.imag
        print("\nS11 = %0.2f dB" % (20*np.log10(np.abs(s11))))
        print("\nR = %0.2f Ohms" % r0)
        print("\nX = %0.2f Ohms" % x0)

        # Move reactance to 50 Ohms
        state = np.argmin(np.abs(statesXs+(x0-50)))
        stateS0 = states[state]
        print("\nSeries capacitance = %0.0f pF"%(statesCs[state]))
        # arduino.write(stateS0+'0000010000')
        s11 = np.array([float(value) for value in vna.readValues("data 0")[idf].split(" ")])
        s11 = s11[0] + s11[1] * 1j
        impedance = 50 * (1. + s11) / (1. - s11)
        r0 = impedance.real
        x0 = impedance.imag
        print("\nR = %0.2f Ohms" % r0)
        print("\nX = %0.2f Ohms" % x0)

        # Move resistance to 50 Ohms
        r0 = 10
        x0 = 50
        a = (r0 - 50)
        b = 2 * 50 * x0
        c = -50 * (r0**2 + x0**2)
        xt0 = (-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
        ct0 = 1 / (2 * np.pi * hw.larmorFreq * 1e6 * xt0)
        state = np.argmin(np.abs(statesCt - ct0))
        stateT0 = states[state]
        # arduino.write(stateS0+stateT0+'10000')
        s11 = np.array([float(value) for value in vna.readValues("data 0")[idf].split(" ")])
        s11 = s11[0] + s11[1] * 1j
        impedance = 50 * (1. + s11) / (1. - s11)
        r0 = impedance.real
        x0 = impedance.imag
        print("\nR = %0.2f Ohms" % r0)
        print("\nX = %0.2f Ohms" % x0)

        # Move reactance to 0
        state = np.argmin(np.abs(statesXm - x0))
        stateM0 = states[state]
        # arduino.write(stateS0+stateT0+stateM0)
        s11 = np.array([float(value) for value in vna.readValues("data 0")[idf].split(" ")])
        s11 = s11[0] + s11[1] * 1j
        impedance = 50 * (1. + s11) / (1. - s11)
        r0 = impedance.real
        x0 = impedance.imag
        print("\nR = %0.2f Ohms" % r0)
        print("\nX = %0.2f Ohms" % x0)

        interface.close()


if __name__ == '__main__':
    seq = AutoTuning()
    seq.sequenceRun()
    seq.sequenceAnalysis(obj='Standalone')
