import configs.hw_config as hw
import seq.mriBlankSeq as blankSeq
import numpy as np
import controller.experiment_gui as ex


class RFOPT(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(RFOPT, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='功率校准', val='rfopt')
        self.addParameter(
            key='larmorFreq', string='中心频率 (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(
            key='rfExAmp', string='RF excitation amplitude (a.u.)', val=0.07, field='RF')
        self.addParameter(
            key='rfExTime', string='激发脉冲宽度 (ms)', val=50.0, field='RF')
        self.addParameter(
            key='rfLobes', string='sinc波峰个数 (0为方波)', val=0, field='RF')
        self.addParameter(key='tau', string='回波间隔 (ms)',
                          val=1, field='SEQ')
        self.addParameter(key='tau2', string='回波间隔2 (ms)', val=1, field='SEQ')
        self.addParameter(key='acqTime', string='采样时长 (ms)',
                          val=2, field='SEQ')
        self.addParameter(key='shimming', string='Shimming',
                          val=[0, 0, 666], field='OTH')
        

    def sequenceInfo(self):
        print("三个90°脉冲时两次回波幅值相同")

    def sequenceTime(self):
        return (0)

    def sequenceRun(self, plotSeq=0):
        larmorFreq = self.mapVals['larmorFreq']  # MHz
        rfExAmp = self.mapVals['rfExAmp']
        rfExTime = self.mapVals['rfExTime']  # us
        shimming = np.array(self.mapVals['shimming'])*1e-4
        shimmingTime = 2e3  # us
        tau = self.mapVals['tau']*1e3  # us
        tau2 = self.mapVals['tau2']*1e3  # us
        nPoints = 256
        acqTime = self.mapVals['acqTime']*1000  # us
        tau_extend = tau - (acqTime + rfExTime)/2
        if tau_extend <= 0:
            print('回波时长不能小于采样时长！')
            return 0
        bw = nPoints / acqTime  # MHz
        nLobes = self.mapVals['rfLobes']
        if (nLobes + 1) % 2 != 0:
            print('Number of lobes must be odd!')
            return 0

        # Initialize the experiment
        samplingPeriod = 1 / bw
        self.expt = ex.Experiment(lo_freq=larmorFreq, rx_t=samplingPeriod)
        samplingPeriod = self.expt.getSamplingRate()
        bw = 1 / samplingPeriod
        acqTime = nPoints / bw
        self.mapVals['acqTime'] = acqTime*1e-3
        self.mapVals['bw'] = bw

        def rfPulse(tStart):
            if nLobes == 0:
                self.rfRecPulse(tStart, rfExTime, rfExAmp)
            else:
                self.rfSincPulse(tStart, rfExTime, rfExAmp, 0, nLobes)

        # Create the sequence
        tim = 20
        self.iniSequence(tim, shimming)
        tim += shimmingTime
        # 第一个脉冲
        rfPulse(tim)
        tim += tau
        # 第二个脉冲
        rfPulse(tim)
        tim += hw.blkTime + rfExTime + tau_extend
        # 第一个回波
        self.rxGateSync(tim, acqTime)
        tim += acqTime + tau2
        # 第三个脉冲
        rfPulse(tim)
        tim += hw.blkTime + rfExTime + tau_extend
        # 第二个回波
        self.rxGateSync(tim, acqTime)
        self.endSequence(1e6)

        if not self.floDict2Exp():
            return 0

        if not plotSeq:
            rxd, msg = self.expt.run()
            print(msg)
            dataFull = self.decimate(rxd['rx0'], 1)
            self.mapVals['data'] = dataFull
            # 将dataFull分成两半，分别对应两个回波
            self.echo1 = dataFull[:dataFull.size]
            self.echo2 = dataFull[dataFull.size:]

        self.expt.__del__()

    def sequenceAnalysis(self):
        tVector = np.linspace(0, self.mapVals['acqTime'], 256)
        self.saveRawData()

        result1 = {
            'widget': 'curve',
            'xData': tVector,
            'yData': [np.abs(self.echo1), np.real(self.echo1), np.imag(self.echo1)],
            'xLabel': 'Time (ms)',
            'yLabel': 'Signal amplitude (mV)',
            'title': 'Signal vs time',
            'legend': ['abs', 'real', 'imag'],
            'row': 0,
            'col': 0
        }
        result2 = {
            'widget': 'curve',
            'xData': tVector,
            'yData': [np.abs(self.echo2), np.real(self.echo2), np.imag(self.echo2)],
            'xLabel': 'Time (ms)',
            'yLabel': 'Signal amplitude (mV)',
            'title': 'Signal vs time',
            'legend': ['abs', 'real', 'imag'],
            'row': 0,
            'col': 1
        }
        return [result1, result2]
