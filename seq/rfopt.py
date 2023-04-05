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
            key='rfAmpMin', string='激发功率下限 (a.u.)', val=0.05, field='RF')
        self.addParameter(
            key='rfAmpMax', string='激发功率上限 (a.u.)', val=0.08, field='RF')
        self.addParameter(key='rfAmpStep', string='激发功率步进 (a.u.)',
                            val=0.001, field='RF')
        self.addParameter(
            key='rfExTime', string='激发脉冲宽度 (μs)', val=50.0, field='RF')
        self.addParameter(
            key='rfLobes', string='sinc波峰个数 (0为方波)', val=0, field='RF')
        self.addParameter(key='tau', string='回波间隔 (ms)',
                          val=1., field='SEQ')
        self.addParameter(key='tau2', string='STE间隔 (ms)', val=1., field='SEQ')
        self.addParameter(key='acqTime', string='采样时长 (ms)',
                          val=0.2, field='SEQ')
        self.addParameter(key='repetitionTime',
                          string='TR(ms)', val=1500, field='SEQ')
        self.addParameter(key='shimming', string='Shimming',
                          val=[300, 300, 1000], field='OTH')
        

    def sequenceInfo(self):
        print("三个90°脉冲时两次回波幅值相同")

    def sequenceTime(self):
        return self.mapVals['repetitionTime']/1000  # sec

    def sequenceRun(self, plotSeq=0):
        larmorFreq = self.mapVals['larmorFreq']  # MHz
        rfExTime = self.mapVals['rfExTime']  # us
        rfAmpMin = self.mapVals['rfAmpMin']
        rfAmpMax = self.mapVals['rfAmpMax']
        rfAmpStep = self.mapVals['rfAmpStep']
        try:
            self.rfExAmp += rfAmpStep
            if self.rfExAmp > rfAmpMax:
                self.rfExAmp = rfAmpMin
        except:
            self.rfExAmp = rfAmpMin
        print('激发功率：', self.rfExAmp)
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
        if nLobes > 0 and nLobes % 2 == 0:
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
                self.rfRecPulse(tStart, rfExTime, self.rfExAmp)
            else:
                self.rfSincPulse(tStart, rfExTime, self.rfExAmp, 0, nLobes)

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
        tim += tau_extend
        # 第二个回波
        self.rxGateSync(tim, acqTime)
        self.endSequence(self.mapVals['repetitionTime']*1000)

        if not self.floDict2Exp():
            return 0

        if not plotSeq:
            rxd, msg = self.expt.run()
            print(msg)
            dataFull = self.decimate(rxd['rx0'], 1)
            self.mapVals['data'] = dataFull
            # 将dataFull分成两半，分别对应两个回波
            self.echo1 = dataFull[:dataFull.size//2]
            self.echo2 = dataFull[dataFull.size//2:]

        self.expt.__del__()

    def sequenceAnalysis(self):
        tVector = np.linspace(0, self.mapVals['acqTime'], self.echo1.size)
        hahnEcho = np.abs(self.echo1)
        ste = np.abs(self.echo2)
        print("回波峰值：", np.max(hahnEcho), np.max(ste))
        self.saveRawData()

        result1 = {
            'widget': 'curve',
            'xData': tVector,
            'yData': [hahnEcho, np.real(self.echo1), np.imag(self.echo1)],
            'xLabel': 'Time (ms)',
            'yLabel': 'Signal amplitude (mV)',
            'title': 'Spin Echo',
            'legend': ['abs', 'real', 'imag'],
            'row': 0,
            'col': 0
        }
        result2 = {
            'widget': 'curve',
            'xData': tVector,
            'yData': [ste, np.real(self.echo2), np.imag(self.echo2)],
            'xLabel': 'Time (ms)',
            'yLabel': 'Signal amplitude (mV)',
            'title': 'Stimulated Echo',
            'legend': ['abs', 'real', 'imag'],
            'row': 0,
            'col': 1
        }
        return [result1, result2]
