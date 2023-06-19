import controller.experiment_gui as ex
import numpy as np
import seq.mriBlankSeq as blankSeq
import configs.hw_config as hw


class RARE(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(RARE, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='3D RARE', val='RARE')
        
        self.addParameter(key='larmorFreq', string='中心频率 (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(key='rfExAmp', string='90°功率', val=0.11, field='RF')
        self.addParameter(key='rfReAmp', string='180°功率', val=0.23, field='RF')
        self.addParameter(key='rfExTime', string='RF excitation time (us)', val=50.0, field='RF')
        self.addParameter(key='rfReTime', string='RF refocus time (us)', val=50.0, field='RF')    
                
        self.addParameter(key='nScans', string='平均次数', val=1, field='IM')
        self.addParameter(key='nPoints', string='像素点数', val=128, field='IM')
        self.addParameter(key='nSlices', string='选层方向点数', val=4, field='IM')
        self.addParameter(key='sliceAmp', string='选层编码步进 (o.u.)', val=0.001, field='IM')
        self.addParameter(key='phaseAmp', string='相位编码步进 (o.u.)', val=0.001, field='IM')
        self.addParameter(key='readAmp', string='读出梯度幅值 (o.u.)', val=0.1, field='IM')

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
        print("============ 3D RARE成像 ============")
        print("我们通过相位编码步进值和持续时间来控制相位编码方向FOV，通过读出梯度值和采样间隔时间(ReadoutTime/nPoints)来控制频率编码方向FOV")

    def sequenceTime(self):  # minutes
        self.sequenceAtributes()
        return self.repetitionTime * self.nShots * self.nSlices * self.nScans / 6e4
        
    def sequenceAtributes(self):
        self.error = True
        super().sequenceAtributes()  # 把mapVals里的键值对读进self里

        read_points = self.nScans*self.nSlices*self.nPoints*self.nPoints*hw.oversamplingFactor
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

        acqTime = self.readoutTime - 2*self.readPadding  # 读出边距
        self.samplingPeriod = acqTime / self.nPoints
        self.nShots = int(self.nPoints / self.etl)
        if self.nPoints % self.etl != 0:
            print('回波链未对齐！')
            return 0

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


    def sequenceRun(self, plotSeq=0, demo=False):
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

        def shot(t_start, segment_idx, slice_amp):
            self.rfRecPulse(t_start, self.rfExTime, self.rfExAmp, 0)  # 激发
            # 读出方向predephase
            t_start += hw.blkTime + self.rfExTime
            gradient(t_start, self.phaseTime, self.ROpreAmp, self.axes[0]) 
        
            for echo_idx in range(self.etl):
                phase_amp = (self.nPoints/2 - (self.etl*segment_idx + echo_idx)) * self.phaseAmp
                tim = t_start + self.rfExTime/2 + self.t_e*(echo_idx + 1) 

                # Refocusing pulse
                self.rfRecPulse(tim - self.t_e/2 - self.rfReTime/2, 
                                self.rfReTime, self.rfReAmp, np.pi/2)
                
                # Slice and Phase encoding
                tim -= self.readoutTime/2 + 3*self.riseTime + self.phaseTime
                gradient(tim, self.phaseTime, phase_amp, self.axes[1])
                gradient(tim, self.phaseTime, slice_amp, self.axes[2])
                
                # Readout
                tim += self.phaseTime + 2 * self.riseTime
                gradient(tim, self.readoutTime, self.readAmp, self.axes[0])

                tim += self.riseTime + self.readPadding
                self.rxGateSync(tim, acq_time)

                tim += acq_time - self.readPadding + self.riseTime
                # phase and slice rewind, no spoil
                gradient(tim, self.phaseTime, -phase_amp, self.axes[1])
                gradient(tim, self.phaseTime, -slice_amp, self.axes[2])

        # --------------------- ↓序列开始↓ ---------------------
        tim = 20
        self.iniSequence(tim, self.shimming)

        for i in range(self.nScans):
            for j in range(self.nSlices):
                slice_amp = (self.nSlices/2 - j) * self.sliceAmp
                for k in range(self.nShots):  # 每次在相位编码方向上采ETL行
                    t0 = 1e5 + ((i*self.nSlices + j)*self.nShots + k)*self.t_r
                    shot(t_start=t0, segment_idx=k, slice_amp=slice_amp)

        self.endSequence(self.nScans * self.nSlices * self.nShots * self.t_r + 2e6)
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
        return True

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

        self.output = self.out = [{
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
        self.saveRawData()
        return self.output
