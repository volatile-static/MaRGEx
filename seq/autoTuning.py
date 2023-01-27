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
import copy
import configs.hw_config as hw

class AutoTuning(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(AutoTuning, self).__init__()
        # Input the parameters
        self.seriesTarget = None
        self.state0 = None
        self.frequencies = None
        self.expt = None
        self.threadpool = None

        # Open arduino serial port
        self.arduino = serial.Serial(port=hw.arduinoPort, baudrate=115200, timeout=.1)
        print('\nArduino connected!')
        time.sleep(1)

        # Initial state
        self.arduino.write('0000000000000000'.encode())

        # Read arduino state
        while self.arduino.in_waiting == 0:
            time.sleep(0.1)
        result = self.arduino.readline()

        # Parameters
        self.addParameter(key='seqName', string='AutoTuningInfo', val='AutoTuning')
        self.addParameter(key='seriesTarget', string='Series target (Ohms)', val=50.0, field='RF')
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

    def runAutoTuning(self):
        start = time.time()
        nCap = 5

        # Combinations
        self.states = ['']*2**nCap
        zero = "00000"
        for state in range(2**nCap):
            prov = bin(state)[2::]
            x = len(prov)
            if x < nCap:
                self.states[state] = zero[0:nCap-x]+prov
            else:
                self.states[state] = prov

        # Series reactance
        cs = np.array([np.Inf, 8, 3.9, 1.8, 1])*1e-9
        self.statesCs = np.zeros(2**nCap)
        self.statesXs = np.zeros(2**nCap)
        for state in range(2**nCap):
            for c in range(nCap):
                if c == 0 and int(self.states[state][0]) == 0:
                    self.statesCs[state] += 0
                else:
                    self.statesCs[state] += int(self.states[state][c])*cs[c]
            if self.statesCs[state]==0:
                self.statesXs[state] = np.Inf
            else:
                self.statesXs[state] = -1/(2*np.pi*hw.larmorFreq*1e6*self.statesCs[state])

        # Tuning capacitor
        ct = np.array([326, 174, 87, 44, 26])*1e-12
        self.statesCt = np.zeros(2**nCap)
        for state in range(2**nCap):
            for c in range(nCap):
                self.statesCt[state] += int(self.states[state][c])*ct[c]

        # Matching capacitors
        cm = np.array([np.Inf, 500, 262, 142, 75])*1e-12
        self.statesCm = np.zeros(2 ** nCap)
        self.statesXm = np.zeros(2 ** nCap)
        for state in range(2 ** nCap):
            for c in range(nCap):
                if c == 0 and int(self.states[state][0]) == 0:
                    self.statesCm[state] += 0
                else:
                    self.statesCm[state] += int(self.states[state][c]) * cm[c]
            if self.statesCm[state]==0:
                self.statesXm[state] = np.Inf
            else:
                self.statesXm[state] = -1/(2*np.pi*hw.larmorFreq*1e6*self.statesCm[state])

        # Open nanoVNA
        # scan serial ports and connect
        interface = get_interfaces()[0]
        interface.open()
        interface.timeout = 0.05
        time.sleep(0.1)
        self.vna = get_VNA(interface)

        # Get frequencies
        self.frequencies = np.array(self.vna.readFrequencies()) * 1e-6
        self.idf = np.argmin(abs(self.frequencies - hw.larmorFreq))

        # Get initial impedance
        state = '10000'
        self.arduino.write((state + '00000100001').encode())
        while self.arduino.in_waiting == 0:
            time.sleep(0.1)
        result = self.arduino.readline()
        s11 = np.array(
            [float(value) for value in self.vna.readValues("data 0")[self.idf].split(" ")])  # "data 0"->S11, "data 1"->S21
        s11 = s11[0] + s11[1] * 1j
        impedance = 50 * (1. + s11) / (1. - s11)
        r0 = impedance.real
        x0 = impedance.imag
        print("\nInput impedance:")
        print("R = %0.2f Ohms" % r0)
        print("X = %0.2f Ohms" % x0)

        if x0 > self.seriesTarget:
            stateCs = self.getCs(0, 17, "1")
        else:
            stateCs = "10000"

        stateCt = self.getCt(0, stateCs, 17, "1")

        stateCm = self.getCm(17, stateCs, stateCt, "1")

        stateCt = self.getCt(stateCt, stateCs, stateCm, "1")

        stateCm = self.getCm(stateCm, stateCs, stateCt, "0")

        print("\nFinal state")
        print(self.states[stateCs]+self.states[stateCt]+self.states[stateCm])

        interface.close()

    def getCs(self, stateCt, stateCm, auto):
        # Sweep series impedances until reactance goes higher than 50 Ohms
        print("\nObtaining series capacitor...")
        n = 0
        x0 = [0.0]
        while x0[-1] < self.seriesTarget and n < 31:
            n += 1
            self.arduino.write((self.states[n] + self.states[stateCt] + self.states[stateCm] + "1").encode())
            while self.arduino.in_waiting == 0:
                time.sleep(0.1)
            result = self.arduino.readline()
            s11 = np.array(
                [float(value) for value in self.vna.readValues("data 0")[self.idf].split(" ")])  # "data 0"->S11, "data 1"->S21
            s11 = s11[0] + s11[1] * 1j
            impedance = 50 * (1. + s11) / (1. - s11)
            x0.append(impedance.imag)

        # Select the value with reactance closest to 50 Ohms
        stateCs = np.argmin(np.abs(np.array(x0) - self.seriesTarget))
        self.arduino.write((self.states[stateCs] + self.states[stateCt] + self.states[stateCm] + auto).encode())
        while self.arduino.in_waiting == 0:
            time.sleep(0.1)
        result = self.arduino.readline()
        s11 = np.array(
            [float(value) for value in self.vna.readValues("data 0")[self.idf].split(" ")])  # "data 0"->S11, "data 1"->S21
        s11 = s11[0] + s11[1] * 1j
        impedance = 50 * (1. + s11) / (1. - s11)
        r0 = impedance.real
        x0 = impedance.imag
        print("Best state:")
        print(self.states[stateCs])
        print("%0.0f pF" % (self.statesCs[stateCs] * 1e12))
        print("S11 = %0.2f dB" % (20 * np.log10(np.abs(s11))))
        print("R = %0.2f Ohms" % r0)
        print("X = %0.2f Ohms" % x0)
        
        return stateCs

    def getCt(self, n0, stateCs, stateCm, auto):
        # Sweep tuning capacitances until resistance goes higher than 50 Ohms
        print("\nObtaining tuning capacitor...")
        n = copy.copy(n0)
        r0 = [0.0]
        while r0[-1] < 50.0 and n < 31:
            n += 1
            self.arduino.write((self.states[stateCs] + self.states[n] + self.states[stateCm] + "1").encode())
            while self.arduino.in_waiting == 0:
                time.sleep(0.1)
            result = self.arduino.readline()
            s11 = np.array(
                [float(value) for value in self.vna.readValues("data 0")[self.idf].split(" ")])  # "data 0"->S11, "data 1"->S21
            s11 = s11[0] + s11[1] * 1j
            impedance = 50 * (1. + s11) / (1. - s11)
            r0.append(impedance.real)

        # Select the value with reactance closest to 50 Ohms
        stateCt = n0 + np.argmin(np.abs(np.array(r0) - 50.0))
        self.arduino.write((self.states[stateCs] + self.states[stateCt] + self.states[stateCm] + auto).encode())
        while self.arduino.in_waiting == 0:
            time.sleep(0.1)
        result = self.arduino.readline()
        s11 = np.array(
            [float(value) for value in self.vna.readValues("data 0")[self.idf].split(" ")])  # "data 0"->S11, "data 1"->S21
        s11 = s11[0] + s11[1] * 1j
        impedance = 50 * (1. + s11) / (1. - s11)
        r0 = impedance.real
        x0 = impedance.imag
        print("Best state:")
        print(self.states[stateCt])
        print("%0.0f pF" % (self.statesCt[stateCt] * 1e12))
        print("S11 = %0.2f dB" % (20*np.log10(np.abs(s11))))
        print("R = %0.2f Ohms" % r0)
        print("X = %0.2f Ohms" % x0)

        return stateCt

    def getCm(self, n0, stateCs, stateCt, auto):
        # Sweep matching capacitances until reactance goes negative
        print("\nObtaining matching capacitor...")
        n = copy.copy(n0)
        x0 = [10000.0]
        while x0[-1] > 0.0 and n > 0:
            n -= 1
            self.arduino.write((self.states[stateCs] + self.states[stateCt] + self.states[n] + "1").encode())
            while self.arduino.in_waiting == 0:
                time.sleep(0.1)
            result = self.arduino.readline()
            s11 = np.array(
                [float(value) for value in self.vna.readValues("data 0")[self.idf].split(" ")])  # "data 0"->S11, "data 1"->S21
            s11 = s11[0] + s11[1] * 1j
            impedance = 50 * (1. + s11) / (1. - s11)
            x0.append(impedance.imag)

        # Select the value with reactance closest to 50 Ohms
        stateCm = n0 - np.argmin(np.abs(np.array(x0)))
        self.arduino.write((self.states[stateCs] + self.states[stateCt] + self.states[stateCm] + auto).encode())
        while self.arduino.in_waiting == 0:
            time.sleep(0.1)
        result = self.arduino.readline()
        s11 = np.array(
            [float(value) for value in self.vna.readValues("data 0")[self.idf].split(" ")])  # "data 0"->S11, "data 1"->S21
        s11 = s11[0] + s11[1] * 1j
        impedance = 50 * (1. + s11) / (1. - s11)
        r0 = impedance.real
        x0 = impedance.imag
        print("Best state:")
        print(self.states[stateCm])
        print("%0.0f pF" % (self.statesCm[stateCm] * 1e12))
        print("S11 = %0.2f dB" % (20*np.log10(np.abs(s11))))
        print("R = %0.2f Ohms" % r0)
        print("X = %0.2f Ohms" % x0)

        return stateCm

if __name__ == '__main__':
    seq = AutoTuning()
    seq.sequenceRun()
    seq.sequenceAnalysis(obj='Standalone')
