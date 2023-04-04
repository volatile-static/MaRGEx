import configs.hw_config as hw
import seq.mriBlankSeq as blankSeq
import numpy as np
import controller.experiment_gui as ex


class GradientEcho(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(GradientEcho, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='梯度回波', val='GradientEcho')
        self.addParameter(
            key='larmorFreq', string='Larmor frequency (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(
            key='rfExAmp', string='RF excitation amplitude (a.u.)', val=0.07, field='RF')
        self.addParameter(
            key='rfExTime', string='RF excitation time (us)', val=50.0, field='RF')
        self.addParameter(key='shimming', string='线性匀场',
                          val=[0, 0, 666], field='OTH')
        # self.addParameter(key='echoTime', string='TE', val=0, field='SEQ')
        self.addParameter(key='gradientChannel', string='梯度通道(x|y|z)', val='x', field='SEQ')
        self.addParameter(key='gradientAmplitude', string='梯度幅度', val=0.1, field='SEQ')
        self.addParameter(key='slewRate', string='梯度斜率', val=0.1, field='SEQ')
        self.addParameter(key='dephaseTime', string='失相位时间', val=0.1, field='SEQ')
        self.addParameter(key='refocusTime', string='重聚时间', val=0.1, field='SEQ')

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
        # self.rfRecPulse(shimmingTime, rfExTime, rfExAmp)
        self.rfSincPulse(shimmingTime, rfExTime, rfExAmp)
        t0 = shimmingTime + hw.blkTime + rfExTime + deadTime
        self.rxGateSync(t0, acqTime)
        self.endSequence(1e6)

        if not self.floDict2Exp():
            return 0

        if not plotSeq:
            rxd, msg = self.expt.run()
            print(msg)
            dataFull = self.decimate(rxd['rx0'], 1)
            self.mapVals['data'] = dataFull

        self.expt.__del__()

    def sequenceAnalysis(self):
        signal = self.mapVals['data']
        bw = self.mapVals['bw']*1e3  # kHz
        nPoints = 100
        deadTime = hw.deadTime*1e-3  # ms
        rfExTime = self.mapVals['rfExTime']*1e-3  # ms
        tVector = np.linspace(rfExTime/2 + deadTime + 0.5/bw,
                              rfExTime/2 + deadTime + (nPoints-0.5)/bw, nPoints)
        fVector = np.linspace(-bw/2, bw/2, nPoints)
        spectrum = np.abs(np.fft.ifftshift(
            np.fft.ifftn(np.fft.ifftshift(signal))))
        fitedLarmor = self.mapVals['larmorFreq'] + \
            fVector[np.argmax(np.abs(spectrum))] * 1e-3
        print('Larmor frequency: %1.5f MHz' % fitedLarmor)
        self.mapVals['signalVStime'] = [tVector, signal]
        self.mapVals['spectrum'] = [fVector, spectrum]
        self.saveRawData()

        # Add time signal to the layout
        result1 = {
            'widget': 'curve',
            'xData': tVector,
            'yData': [np.abs(signal), np.real(signal), np.imag(signal)],
            'xLabel': 'Time (ms)',
            'yLabel': 'Signal amplitude (mV)',
            'title': 'Signal vs time',
            'legend': ['abs', 'real', 'imag'],
            'row': 0,
            'col': 0
        }

        # Add frequency spectrum to the layout
        result2 = {
            'widget': 'curve',
            'xData': fVector,
            'yData': [spectrum],
            'xLabel': 'Frequency (kHz)',
            'yLabel': 'Spectrum amplitude (a.u.)',
            'title': 'Spectrum',
            'legend': [''],
            'row': 1,
            'col': 0
        }
        return [result1, result2]
