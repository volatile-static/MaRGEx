from cmath import phase
from re import S
import controller.experiment_gui as ex
import numpy as np
import seq.mriBlankSeq as blankSeq
import configs.hw_config as hw


class A2RE(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(A2RE, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='A2RE', val='A2RE')
        
        self.addParameter(key='larmorFreq', string='中心频率 (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(key='rfExAmp', string='90°功率', val=0.11, field='RF')
        self.addParameter(key='rfReAmp', string='180°功率', val=0.23, field='RF')
        self.addParameter(key='rfExTime', string='RF excitation time (us)', val=50.0, field='RF')
        self.addParameter(key='rfReTime', string='RF refocus time (us)', val=50.0, field='RF')    
                
        self.addParameter(key='nPoints', string='像素点数', val=128, field='IM')
        self.addParameter(key='nSlices', string='选层方向点数', val=4, field='IM')
        self.addParameter(key='sliceAmp', string='选层编码幅值 (o.u.)', val=0.01, field='IM')
        self.addParameter(key='phaseAmp', string='相位编码幅值 (o.u.)', val=0.1, field='IM')
        self.addParameter(key='readAmp', string='频率编码幅值 (o.u.)', val=0.1, field='IM')

        self.addParameter(key='etl', string='回波链长度', val=4, field='SEQ')
        self.addParameter(key='echoSpacing', string='回波间隔 (ms)', val=5, field='SEQ')
        self.addParameter(key='repetitionTime',string='TR (ms)', val=2000.0, field='SEQ')
        self.addParameter(key='phaseTime', string='相位编码时长 (ms)', val=2, field='SEQ')
        self.addParameter(key='readoutTime', string='读出时长 (ms)', val=3.0, field='SEQ')
        self.addParameter(key='readPadding', string='读出边距 (μs)', val=1.0, field='SEQ')

        self.addParameter(key='shimming', string='线性匀场 [x,y,z]', val=[210.0, 210.0, 895.0], field='OTH')
        self.addParameter(key='axes', string='[读出，相位，选层]', val=[0, 1, 2], field='OTH')
        self.addParameter(key='preEmphasisFactor', string='预加重比例', val=1.0, field='OTH')

    def sequenceInfo(self):
        print("============ Averaging Acquisition with Relaxation Enhancement ============")
        print("利用RARE 3D的回波链进行平均的SE序列")

    def sequenceTime(self):  # minutes
        self.sequenceAtributes()
        return self.repetitionTime * self.nPoints * self.nSlices / 6e4
        
    def sequenceAtributes(self):
        self.error = True
        super().sequenceAtributes()  # 把mapVals里的键值对读进self里

        read_points = self.etl*self.nSlices*self.nPoints*self.nPoints*hw.oversamplingFactor
        print('采样深度: ', read_points)
        if read_points > hw.maxRdPoints:
            print('读出点数过多！')
            # return 0
        
        self.riseTime = 300
        self.riseSteps = 20
            
        self.shimming = np.array(self.shimming) / 1e4
        self.phaseTime *= 1e3  # μs
        self.readoutTime *= 1e3  # μs
        self.t_r = self.repetitionTime * 1e3  # μs
        self.t_e = self.echoSpacing * 1e3  # μs
        self.tau = self.t_e / 2  # μs

        acqTime = self.readoutTime - 2*self.readPadding  # 读出边距
        self.samplingPeriod = acqTime / self.nPoints
        self.phaseGrads = np.linspace(-self.phaseAmp, self.phaseAmp, self.nPoints, False)  # 相位编码梯度
        self.sliceGrads = np.linspace(-self.sliceAmp, self.sliceAmp, self.nSlices, False)  # 选层编码梯度

        echoSpacingMin = self.rfReTime + hw.blkTime + 2*self.phaseTime + 6*self.riseTime + self.readoutTime
        print('最小回波间隔: ', echoSpacingMin, 'μs')
        if self.t_e < echoSpacingMin:
            self.t_e = echoSpacingMin + 2
            self.mapVals['echoSpacing'] = self.t_e / 1e3
            
        # 计算ReadOut predephase梯度大小
        self.ROpreAmp = self.preEmphasisFactor * 0.5 * self.readAmp * (self.readoutTime + self.riseTime) / (self.phaseTime + self.riseTime)
        if np.abs(self.ROpreAmp) > 1:
            print('ReadOut predephase梯度过大！')
            return 0

        self.phaseAmpMax = self.phaseAmp * (self.nPoints - 1) / 2
        if self.phaseAmpMax > 1:
            print('相位编码过大！')
            return 0
        self.error = False


    def sequenceRun(self, plotSeq=0):
        if self.error:
            return 0
        self.error = True

        self.expt = ex.Experiment(lo_freq=self.mapVals['larmorFreq'], rx_t=self.samplingPeriod)
        self.mapVals['samplingRate'] = self.expt.getSamplingRate()  # 采样间隔
        acq_time = self.mapVals['samplingRate'] * self.nPoints
        print('采样率：', 1e3/self.mapVals['samplingRate'], ' (kHz)')
        
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

        def shot(t_start, slice_amp, phase_amp):
            self.rfRecPulse(t_start, self.rfExTime, self.rfExAmp)  # 激发

            # 读出方向predephase
            t_start += hw.blkTime + self.rfExTime
            gradient(t_start, self.phaseTime, self.ROpreAmp, self.axes[0]) 
        
            t_start += self.tau - self.rfExTime/2
            for i in range(self.etl):
                t_refocus = t_start + self.t_e*i 
                t_echo = t_refocus + self.tau

                # Refocusing pulse
                self.rfRecPulse(t_refocus - hw.blkTime - self.rfReTime/2, 
                                self.rfReTime, self.rfReAmp, np.pi/2)
                
                # Readout
                t_read = t_echo - self.readoutTime/2
                gradient(t_read - self.riseTime, self.readoutTime, self.readAmp, self.axes[0])
                self.rxGateSync(t_read + self.readPadding, acq_time)

                # Slice and Phase encoding
                t_phase = t_read - 3*self.riseTime - self.phaseTime
                gradient(t_phase, self.phaseTime, phase_amp, self.axes[1])
                gradient(t_phase, self.phaseTime, slice_amp, self.axes[2])
                
                # phase and slice rewind, no spoil
                t_rewind = t_read + self.readoutTime + self.riseTime + 1
                gradient(t_rewind, self.phaseTime, -phase_amp, self.axes[1])
                gradient(t_rewind, self.phaseTime, -slice_amp, self.axes[2])

        # --------------------- ↓序列开始↓ ---------------------
        tim = 20
        self.iniSequence(tim, self.shimming)

        for i in range(self.nSlices):
            for j in range(self.nPoints): 
                shot(
                    t_start=1e5 + (i*self.nPoints + j)*self.t_r, 
                    slice_amp=self.sliceGrads[i], 
                    phase_amp=self.phaseGrads[j]
                )

        self.endSequence(self.nSlices * self.nShots * self.t_r + 2e6)
        # --------------------- ↑序列结束↑ ---------------------

        if not self.floDict2Exp():  # 验证时序
            return 0

        if not plotSeq:
            # Run the experiment and get data
            rxd, msgs = self.expt.run()
            print(msgs)
            self.mapVals['dataOver'] = rxd['rx0']

            self.expt.__del__()
        self.error = False

    def sequenceAnalysis(self):
        if self.error:
            self.saveRawData()
            return []
        
        # 对每次读出分别降采样
        data_full = self.decimate(self.mapVals['dataOver'], self.etl * self.nPoints * self.nSlices)

        # 对回波链取平均
        data_raw = np.mean(np.reshape(data_full, (self.nSlices, self.nPoints, self.etl, -1)), 2)
        
        ksp = self.mapVals['ksp3d'] = np.reshape(data_raw, (self.nSlices, self.nPoints, self.nPoints))
        img = self.mapVals['img3d'] = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(ksp)))  # 重建
        
        self.output = [{
            'widget': 'image',
            'data': np.concatenate((np.abs(img), np.angle(img))),
            'xLabel': '相位编码',
            'yLabel': '频率编码',
            'title': '幅值图与相位图',
            'row': 0,
            'col': 0
        }, {
            'widget': 'image',
            'data': np.log(np.abs(ksp)),
            'xLabel': 'mV',
            'yLabel': 'ms',
            'title': 'k-Space',
            'row': 0,
            'col': 1
        }]
        self.saveRawData()
        return self.output
