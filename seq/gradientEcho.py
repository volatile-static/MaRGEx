from turtle import st
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
                          val=[300, 300, 1000], field='OTH')
        self.addParameter(
            key='slewRate', string='梯度斜率 (us/o.u.)', val=1000, field='OTH')
        self.addParameter(
            key='stepsRate', string='梯度步进速率 (step/o.u.)', val=200, field='OTH')
        self.addParameter(key='gradUpdateRate', string='梯度更新速率 (o.u./s)',
                          val=0.1, field='OTH')
        self.addParameter(key='dephaseTime', string='dephase time (us)',
                          val=100, field='SEQ')
        self.addParameter(key='gradientChannel',
                          string='梯度通道 (0/1/2)', val=0, field='SEQ')
        self.addParameter(key='gradientAmplitude',
                          string='梯度幅度', val=0.1, field='SEQ')

    def sequenceInfo(self):
        print("单纯的梯度回波")

    def sequenceTime(self):
        return (0)

    def sequenceRun(self, plotSeq=0):
        # us/o.u., slew rate for gradient rises
        slewRate = self.mapVals['slewRate']
        # steps/o.u., steps rate for gradient rises
        stepsRate = self.mapVals['stepsRate']
        larmorFreq = self.mapVals['larmorFreq']  # MHz
        rfExAmp = self.mapVals['rfExAmp']
        rfExTime = self.mapVals['rfExTime']  # us
        shimming = np.array(self.mapVals['shimming'])*1e-4
        shimmingTime = 2e5  # us
        dephaseTime = self.mapVals['dephaseTime']  # us
        gradientAmplitude = self.mapVals['gradientAmplitude']
        raiseTime = slewRate * gradientAmplitude
        refocusTime = raiseTime + 2*dephaseTime
        acqTime = refocusTime + raiseTime*2  # us
        gSteps = int(np.round(gradientAmplitude * stepsRate))
        if gSteps < 1:
            print("梯度步进过小！")
            gSteps = 1
        nPoints = 256
        bw = nPoints / acqTime  # MHz

        # Initialize the experiment
        samplingPeriod = 1 / bw
        self.expt = ex.Experiment(
            lo_freq=larmorFreq,
            rx_t=samplingPeriod,
            init_gpa=True,
            grad_max_update_rate=self.mapVals['gradUpdateRate']
        )
        samplingPeriod = self.expt.getSamplingRate()
        bw = 1 / samplingPeriod
        acqTime = nPoints / bw
        self.mapVals['acqTime'] = acqTime*1e-3
        self.mapVals['bw'] = bw

        # Create the sequence
        tim = 20
        self.iniSequence(tim, shimming)
        tim += shimmingTime
        self.rfRecPulse(tim, rfExTime, rfExAmp)
        tim += hw.blkTime + rfExTime + hw.deadTime
        self.gradTrap(tim, raiseTime, dephaseTime, -gradientAmplitude,
                      gSteps, self.mapVals['gradientChannel'], shimming)
        tim += dephaseTime + 2*raiseTime + 1
        self.gradTrap(tim, raiseTime, refocusTime, gradientAmplitude,
                      gSteps, self.mapVals['gradientChannel'], shimming)
        self.rxGateSync(tim, acqTime)
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
        nPoints = 256
        tVector = np.linspace(0, nPoints/bw, nPoints)
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
