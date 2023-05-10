import controller.experiment_gui as ex
import numpy as np
import seq.mriBlankSeq as blankSeq
import configs.hw_config as hw


class T1SE(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(T1SE, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='T1w Spin Echo', val='T1SE')
        self.addParameter(
            key='larmorFreq', string='中心频率 (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(
            key='rfExAmp', string='90°功率', val=0.3, field='RF')
        self.addParameter(
            key='rfReAmp', string='180°功率', val=0.6, field='RF')
        self.addParameter(
            key='rfExTime', string='RF excitation time (us)', val=30.0, field='RF')
        self.addParameter(key='repetitionTime',
                          string='Repetition time (ms)', val=600., field='SEQ')
        self.addParameter(
            key='nPoints', string='Number of points', val=256, field='SEQ')
        self.addParameter(key='phaseAmp', string='相位编码幅值', val=0.1, field='IM')
        self.addParameter(key='readAmp', string='读出梯度幅值', val=0.1, field='IM')
        self.addParameter(key='preparationAmp',
                          string='预读出幅值', val=0.1, field='IM')
        self.addParameter(key='shimming', string='线性匀场',
                          val=[300, 300, 1000], field='OTH')

    def sequenceInfo(self):
        print(" ")

    def sequenceTime(self):  # minutes, scanTime
        return self.mapVals['repetitionTime'] * self.mapVals['nPoints'] / 6e4

    def sequenceRun(self, plotSeq=0):
        # 读取输入的参数
        larmorFreq = self.mapVals['larmorFreq']  # MHz
        rfExAmp = self.mapVals['rfExAmp']
        rfReAmp = self.mapVals['rfReAmp']
        rfExTime = self.mapVals['rfExTime']  # us
        tRepetition = self.mapVals['repetitionTime']*1e3  # us
        nPoints = self.mapVals['nPoints']
        shimming = np.array(self.mapVals['shimming'])*1e-4
        phaseAmp = self.mapVals['phaseAmp']  # 1.6mV
        readAmp = self.mapVals['readAmp']  # 25mV
        preparationAmp = self.mapVals['preparationAmp']  # 70mV

        # 定义时序（μs）
        samplingPeriod = 30
        echoSpacing = 14700
        tShimming = 2000
        tRamp = 300
        tPhase = 2000

        channelPhase = 0
        channelRead = 1

        # 校验序列参数
        if nPoints * phaseAmp > 1:
            print("Phase encoding amplitude is too large!")
            return

        # Initialize the experiment
        self.expt = ex.Experiment(lo_freq=larmorFreq, rx_t=samplingPeriod)
        samplingPeriod = self.expt.getSamplingRate()
        acqTime = nPoints * samplingPeriod

        # --------------------- ↓序列开始↓ ---------------------
        self.iniSequence(20, shimming)
        for i in range(nPoints):
            t0 = tShimming + i*tRepetition
            # Excitation pulse
            self.rfRecPulse(t0, rfExTime, rfExAmp)
            self.rfRecPulse(t0 + echoSpacing/2, rfExTime, rfReAmp)
            # self.rfSincPulse(t0, rfExTime, rfExAmp, 0, 3)
            # self.rfSincPulse(t0 + echoSpacing/2, rfExTime, rfReAmp, 0, 3)

            # phase encoding
            gPhAmp = (nPoints/2 - i) * phaseAmp
            self.gradTrap(
                tStart=t0 + rfExTime + tRamp,
                gRiseTime=tRamp,
                gFlattopTime=tPhase,
                gAmp=gPhAmp,
                gSteps=int(hw.stepsRate * abs(gPhAmp)) + 1,
                gAxis=channelPhase,
                shimming=shimming
            )

            # read out
            self.gradTrap(
                tStart=t0 + rfExTime + tRamp,
                gRiseTime=tRamp,
                gFlattopTime=tPhase,
                gAmp=preparationAmp,
                gSteps=int(hw.stepsRate * preparationAmp),
                gAxis=channelRead,
                shimming=shimming
            )
            t1 = t0 + echoSpacing + (rfExTime - acqTime)/2
            self.rxGateSync(t1, acqTime)

            self.setGradientRamp(
                tStart=t1 - tRamp,
                gradRiseTime=tRamp,
                nStepsGradRise=int(hw.stepsRate * readAmp),
                g0=0,
                gf=readAmp,
                gAxes=channelRead,
                shimming=shimming
            )
            self.setGradientRamp(
                tStart=t1 + acqTime,
                gradRiseTime=tRamp,
                nStepsGradRise=int(hw.stepsRate * abs(preparationAmp - readAmp)),
                g0=readAmp,
                gf=preparationAmp,
                gAxes=channelRead,
                shimming=shimming
            )
            self.setGradientRamp(
                tStart=t1 + acqTime + tRamp + tPhase,
                gradRiseTime=tRamp,
                nStepsGradRise=int(hw.stepsRate * preparationAmp),
                g0=preparationAmp,
                gf=0,
                gAxes=channelRead,
                shimming=shimming
            )

        self.endSequence(nPoints*tRepetition)
        # --------------------- ↑序列结束↑ ---------------------

        if not self.floDict2Exp():  # 验证时序
            return 0

        if not plotSeq:
            # Run the experiment and get data
            rxd, msgs = self.expt.run()

            # Decimate the signal
            dataFull = self.decimate(rxd['rx0'], nPoints)
            dataFull *= hw.adcFactor  # Here I normalize to get the result in mV
            kSpace = np.reshape(dataFull, (1, nPoints, nPoints))
            img = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(kSpace)))  # 重建

            self.mapVals['data'] = kSpace
            self.mapVals['img'] = img

        self.expt.__del__()

    def sequenceAnalysis(self, obj=''):
        self.saveRawData()
        return [
            {
                'widget': 'image',
                'data': np.abs(self.mapVals['img']),
                'xLabel': 'xLabel',
                'yLabel': 'yLabel',
                'title': 'Image',
                'row': 0,
                'col': 0
            },
            {
                'widget': 'image',
                'data': np.log10(np.abs(self.mapVals['data'])),
                'xLabel': 'xLabel',
                'yLabel': 'yLabel',
                'title': 'k-Space',
                'row': 0,
                'col': 1
            }
        ]
