"""
@author: J.M. Algarín, february 03th 2022
MRILAB @ I3M
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
import time
import experiment as ex
import numpy as np
import seq.mriBlankSeq as blankSeq  # Import the mriBlankSequence for any new sequence.
import scipy.signal as sig
import configs.hw_config as hw

class Noise(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(Noise, self).__init__()
        # Input the parameters
        self.rxChannel = None
        self.nPoints = None
        self.bw = None
        self.freqOffset = None
        self.addParameter(key='seqName', string='NoiseInfo', val='Noise')
        self.addParameter(key='freqOffset', string='RF frequency offset (kHz)', val=0.0, field='RF')
        self.addParameter(key='nPoints', string='Number of points', val=2500, field='RF')
        self.addParameter(key='bw', string='Acquisition bandwidth (kHz)', val=50.0, field='RF')
        self.addParameter(key='rxChannel', string='Rx channel', val=0, field='RF')

    def sequenceInfo(self):
        print(" ")
        print("Noise")
        print("Author: Dr. J.M. Algarín")
        print("Contact: josalggui@i3m.upv.es")
        print("mriLab @ i3M, CSIC, Spain")
        print("Get a noise measurement")
        print(" ")

    def sequenceTime(self):
        return(0)  # minutes, scanTime

    def sequenceRun(self, plotSeq=0, demo=False):
        init_gpa = False

        # Create the inputs automatically as class properties
        for key in self.mapKeys:
            setattr(self, key, self.mapVals[key])

        # Fix units to MHz and us
        self.freqOffset *= 1e-3 # MHz
        self.bw *= 1e-3 # MHz

        if demo:
            dataR = np.random.randn(self.nPoints*hw.oversamplingFactor)
            dataC = np.random.randn(self.nPoints*hw.oversamplingFactor)
            data = dataR+1j*dataC
            data = sig.decimate(data, hw.oversamplingFactor, ftype='fir', zero_phase=True)
            acqTime = self.nPoints/self.bw
            tVector = np.linspace(0, acqTime, num=self.nPoints) * 1e-3  # ms
            spectrum = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(data)))
            fVector = np.linspace(-self.bw / 2, self.bw / 2, num=self.nPoints) * 1e3  # kHz
            self.dataTime = [tVector, data]
            self.dataSpec = [fVector, spectrum]
            time.sleep(0.5)
        else:
            self.bw = self.bw * hw.oversamplingFactor
            samplingPeriod = 1 / self.bw
            self.expt = ex.Experiment(lo_freq=hw.larmorFreq + self.freqOffset*1e-3,
                                      rx_t=samplingPeriod,
                                      init_gpa=init_gpa,
                                      gpa_fhdo_offset_time=(1 / 0.2 / 3.1),
                                      print_infos=False)
            samplingPeriod = self.expt.get_rx_ts()[0]
            self.bw = 1/samplingPeriod/hw.oversamplingFactor
            acqTime = self.nPoints/self.bw

            # SEQUENCE
            self.iniSequence(20, np.array((0, 0, 0)))
            self.rxGate(30, acqTime, channel=self.rxChannel)
            self.endSequence(acqTime+40)
            if self.floDict2Exp():
                pass
            else:
                return 0

            if plotSeq == 0:
                rxd, msgs = self.expt.run()
                data = sig.decimate(rxd['rx%i' % self.rxChannel]*hw.adcFactor, hw.oversamplingFactor, ftype='fir', zero_phase=True)
                self.mapVals['data'] = data
                tVector = np.linspace(0, acqTime, num=self.nPoints) * 1e-3  # ms
                spectrum = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(data)))
                fVector = np.linspace(-self.bw / 2, self.bw / 2, num=self.nPoints) * 1e3  # kHz
                self.dataTime = [tVector, data]
                self.dataSpec = [fVector, spectrum]
            self.expt.__del__()

    def sequenceAnalysis(self, obj=''):
        noiserms = np.std(self.dataTime[1])
        self.mapVals['RMS noise'] = noiserms
        self.mapVals['sampledPoint'] = noiserms # for sweep method
        self.saveRawData()
        print('\nrms noise: %0.5f mV' % noiserms)

        # Plot signal versus time
        result1 = {'widget': 'curve',
                   'xData': self.dataTime[0],
                   'yData': [np.abs(self.dataTime[1]), np.real(self.dataTime[1]), np.imag(self.dataTime[1])],
                   'xLabel': 'Time (ms)',
                   'yLabel': 'Signal amplitude (mV)',
                   'title': 'Noise vs time',
                   'legend': ['abs', 'real', 'imag'],
                   'row': 0,
                   'col': 0}

        # Plot spectrum
        result2 = {'widget': 'curve',
                   'xData': self.dataSpec[0],
                   'yData': [np.abs(self.dataSpec[1])],
                   'xLabel': 'Frequency (kHz)',
                   'yLabel': 'Mag FFT (a.u.)',
                   'title': 'Noise spectrum',
                   'legend': [''],
                   'row': 1,
                   'col': 0}

        self.out = [result1, result2]

        return self.out


if __name__=='__main__':
    # seq = Noise()
    # seq.sequenceRun()
    # seq.sequenceAnalysis(obj='Standalone')

    import pyqtgraph.examples
    pyqtgraph.examples.run()

