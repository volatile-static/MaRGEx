import configs.hw_config as hw
import scipy.signal as sig
import seq.mriBlankSeq as blankSeq
import numpy as np
import controller.experiment_gui as ex
import os
import sys

# *****************************************************************************
# Add path to the working directory
path = os.path.realpath(__file__)
ii = 0
for char in path:
    if (char == '\\' or char == '/') and path[ii+1:ii+14] == 'PhysioMRI_GUI':
        sys.path.append(path[0:ii+1]+'PhysioMRI_GUI')
        sys.path.append(path[0:ii+1]+'marcos_client')
    ii += 1
# ******************************************************************************


class SoftFID(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(SoftFID, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='FIDinfo', val='SoftFID')
        self.addParameter(key='larmorFreq', string='Larmor frequency (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(key='rfExAmp', string='RF excitation amplitude (a.u.)', val=0.3, field='RF')
        self.addParameter(key='rfExTime', string='RF excitation time (us)', val=30.0, field='RF')
        self.addParameter(key='shimming', string='Shimming', val=[0, 0, 666], field='OTH')


    def sequenceInfo(self):
        print("用软脉冲的FID")


    def sequenceTime(self):
        return (0)


    def sequenceRun(self, plotSeq=0):
        larmorFreq = self.mapVals['larmorFreq']  # MHz
        rfExAmp = self.mapVals['rfExAmp']
        rfExTime = self.mapVals['rfExTime']  # us
        deadTime = hw.deadTime
        shimming = np.array(self.mapVals['shimming'])*1e-4
        shimmingTime = 2e3  # us
        nPoints = 100
        acqTime = 1e3  # us
        bw = nPoints / acqTime  # MHz

        # Initialize the experiment
        samplingPeriod = 1 / bw
        self.expt = ex.Experiment(lo_freq=larmorFreq, rx_t=samplingPeriod)
        samplingPeriod = self.expt.getSamplingRate()
        bw = 1 / samplingPeriod
        acqTime = nPoints / bw
        self.mapVals['acqTime'] = acqTime*1e-3
        self.mapVals['bw'] = bw

        # Create the sequence
        self.iniSequence(20, shimming)
        self.rfRecPulse(shimmingTime, rfExTime, rfExAmp)
        t0 = shimmingTime + hw.blkTime + rfExTime + deadTime
        self.rxGateSync(t0, acqTime)
        self.endSequence(t0 + acqTime + 1)

        if not self.floDict2Exp():
            return 0
        
        if not plotSeq:
            rxd, msg = self.expt.run()
            print(msg)
            dataFull = self.decimate(rxd['rx0'], 1)
            self.mapVals['data'] = dataFull

        self.expt.__del__()


    def sequenceAnalysis(self):
        return []
