import numpy as np
from numpy.random import uniform

import configs.hw_config as hw
import controller.experiment_gui as ex
import seq.mriBlankSeq as blankSeq


class GRE3D(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(GRE3D, self).__init__()
        # Input the parameters
        self.error = False
        self.addParameter(key='seqName', string='梯度回波成像', val='GRE2D')

        self.addParameter(key='nScans', string='平均次数', val=1, field='IM')
        self.addParameter(key='nPoints', string='像素点数', val=128, field='IM')
        self.addParameter(key='nSlices', string='选层方向点数', val=4, field='IM')
        self.addParameter(key='sliceAmp', string='选层编码步进 (o.u.)', val=0.001, field='IM')
        self.addParameter(key='phaseAmp', string='相位编码步进 (o.u.)', val=0.001, field='IM')
        self.addParameter(key='readAmp', string='读出梯度幅值 (o.u.)', val=0.1, field='IM')

        self.addParameter(key='rfExAmp', string='激发功率 (a.u.)', val=0.05, field='RF')
        self.addParameter(key='rfExTime', string='激发时长 (μs)', val=50.0, field='RF')
        self.addParameter(key='larmorFreq', string='中心频率 (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(key='bwCalib', string='调频带宽 (MHz)', val=0.02, field='RF')

        self.addParameter(key='echoSpacing', string='TE (ms)', val=2.5, field='SEQ')
        self.addParameter(key='repetitionTime', string='TR (ms)', val=100, field='SEQ')
        self.addParameter(key='phaseTime', string='相位编码时长 (ms)', val=0.3, field='SEQ')
        self.addParameter(key='readoutTime', string='读出时长 (ms)', val=1.0, field='SEQ')
        self.addParameter(key='readPadding', string='读出边距 (μs)', val=1.0, field='SEQ')

        self.addParameter(key='shimming', string='线性匀场 [x,y,z]', val=[210.0, 210.0, 895.0], field='OTH')
        self.addParameter(key='axes', string='[读出，相位，选层]', val=[0, 1, 2], field='OTH')
        self.addParameter(key='riseTime', string='梯度上升时间 (μs)', val=300, field='OTH')
        self.addParameter(key='riseSteps', string='梯度上升步数', val=20, field='OTH')

    def sequenceInfo(self):
        print("============ 梯度回波成像 ============")
        print("我们通过相位编码步进值和持续时间来控制相位编码方向FOV，通过读出梯度值和采样间隔时间(ReadoutTime/nPoints)来控制频率编码方向FOV")

    def sequenceTime(self):
        return self.mapVals['repetitionTime'] * self.mapVals['nPoints'] * self.mapVals['nSlices'] * self.mapVals['nScans'] / 6e4

    def sequenceAtributes(self):
        self.error = True
        super().sequenceAtributes()  # 把mapVals里的键值对读进self里

        read_points = self.nScans*self.nSlices*self.nPoints*self.nPoints*hw.oversamplingFactor
        print('采样深度: ', read_points)
        if read_points > hw.maxRdPoints:
            print('读出点数过多！')
            return 0

        self.spoilDelay = 100
        self.shimming = np.array(self.shimming) / 1e4
        self.phaseTime *= 1e3  # μs
        self.readoutTime *= 1e3  # μs
        self.t_e = self.echoSpacing * 1e3  # μs
        self.t_r = self.repetitionTime * 1e3  # μs

        acqTime = self.readoutTime - 2*self.readPadding  # 读出边距
        self.samplingPeriod = acqTime / self.nPoints

        # 计算ReadOut predephase梯度大小
        self.ROpreAmp = -0.5 * self.readAmp * (self.readoutTime + self.riseTime) / (self.phaseTime + self.riseTime)
        if np.abs(self.ROpreAmp) > 1:
            print('ReadOut predephase梯度过大！')
            return 0

        self.phaseAmpMax = self.phaseAmp * (self.nPoints - 1) / 2
        if self.phaseAmpMax > 1:
            print('相位编码过大！')
            return 0
        self.error = False


    def sequenceRun(self, plot_seq=0):
        if self.error:
            return 0
        self.error = True

        if self.bwCalib > 0:  # 自动调频
            hw.larmorFreq = self.freqCalibration(self.bwCalib, 0.001)
            hw.larmorFreq = self.freqCalibration(self.bwCalib)
            self.mapVals['larmorFreq'] = hw.larmorFreq

        def gradient(t, flat, amp, channel):
            self.gradTrap(
                tStart=t,
                gRiseTime=self.riseTime,
                gFlattopTime=flat,
                gAmp=amp,
                gSteps=self.riseSteps,
                gAxis=channel,
                shimming=self.shimming
            )

        self.expt = ex.Experiment(lo_freq=self.mapVals['larmorFreq'], rx_t=self.samplingPeriod)
        self.mapVals['samplingRate'] = self.expt.getSamplingRate()  # 采样间隔
        acq_time = self.mapVals['samplingRate'] * self.nPoints
        print('采样率：', 1e3/self.mapVals['samplingRate'], ' (kHz)')

        # --------------------- ↓序列开始↓ ---------------------
        tim = 20
        self.iniSequence(tim, self.shimming)

        for i in range(self.nScans):
            for j in range(self.nSlices):
                slice_amp = (self.nSlices/2 - j) * self.sliceAmp
                for k in range(self.nPoints):
                    phase_amp = (self.nPoints/2 - k) * self.phaseAmp
                    rf_ex_phase = 0#.5*117*((j+k)*(j+k) + (j+k) + 2)%360

                    tim = 1e5 + (i*self.nPoints*self.nSlices + j*self.nPoints + k) * self.t_r
                    self.rfRecPulse(tim, self.rfExTime, self.rfExAmp, rf_ex_phase / 180 * np.pi)
                    
                    tim += self.rfExTime
                    # 相位编码、选层编码和读出方向predephase同时开
                    gradient(tim, self.phaseTime, phase_amp, self.axes[1])
                    gradient(tim, self.phaseTime, slice_amp, self.axes[2])
                    gradient(tim, self.phaseTime, self.ROpreAmp, self.axes[0]) 
                
                    tim += self.phaseTime + 2 * self.riseTime + 1
                    gradient(tim, self.readoutTime, self.readAmp, self.axes[0])

                    tim += self.riseTime + self.readPadding
                    self.rxGateSync(tim, acq_time)

                    # slice and spoil
                    tim += acq_time - self.readPadding + self.riseTime + self.spoilDelay
                    gradient(tim, self.phaseTime, self.phaseAmpMax*uniform(-1, 1), self.axes[2])
                    # phase rewind
                    gradient(tim, self.phaseTime, -phase_amp, self.axes[1])

        self.endSequence(self.nSlices * self.nPoints * self.nScans * self.t_r + 2e6)
        # --------------------- ↑序列结束↑ ---------------------

        if not self.floDict2Exp():
            return 0

        if not plot_seq:
            rxd, msg = self.expt.run()
            print(msg)
            self.mapVals['dataOver'] = rxd['rx0']

        self.expt.__del__()
        self.error = False

    def sequenceAnalysis(self):
        if self.error:
            self.saveRawData()
            return []
        data_over = np.reshape(self.mapVals['dataOver'], (self.mapVals['nScans'], -1))
        print('数据维度：', data_over.shape)
        data_average = np.average(data_over, axis=0).reshape((self.nSlices, -1))
        data_full = []
        for i in range(self.nSlices):
            data_full.append(self.decimate(data_average[i], self.nPoints))
        ksp = self.mapVals['ksp3d'] = np.reshape(data_full, (self.nSlices, self.nPoints, self.nPoints))
        img = self.mapVals['img3d'] = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(ksp)))  # 重建
        self.saveRawData()

        return [{
            'widget': 'image',
            'data': np.abs(img),
            'xLabel': '相位编码',
            'yLabel': '频率编码',
            'title': '幅值图',
            'row': 0,
            'col': 0
        }, {
            'widget': 'image',
            'data': np.log10(np.abs(ksp)),
            'xLabel': 'mV',
            'yLabel': 'ms',
            'title': 'k-Space',
            'row': 0,
            'col': 1
        }]
