import numpy as np
from numpy.fft import ifftn, ifftshift
from scipy.signal import decimate
from experiment import Experiment
import configs.hw_config as hw
import seq.mriBlankSeq as blankSeq


class A2RE(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(A2RE, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='A2RE', val='A2RE')

        self.addParameter(key='larmorFreq', string='中心频率 (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(key='rfExAmp', string='90°功率', val=0.11, field='RF')
        self.addParameter(key='rfReAmp', string='180°功率', val=0.23, field='RF')
        self.addParameter(key='rfExTime', string='RF excitation time (us)', val=50.0, field='RF')

        self.addParameter(key='repetitionTime', string='TR (ms)', val=2000.0, field='IM')
        self.addParameter(key='axesOrientation', string='方位 [rd, ph, sl]', val=[0, 1, 2], field='IM')
        self.addParameter(key='nPoints', string='点数 [rd, ph, sl]', val=[128, 128, 4], field='IM')
        self.addParameter(key='voxel', string='体素 (mm)', val=[1.0, 1.0, 2.0], field='IM')

        self.addParameter(key='etl', string='回波链长度', val=4, field='SEQ')
        self.addParameter(key='echoSpacing', string='回波间隔 (ms)', val=5, field='SEQ')
        self.addParameter(key='phaseTime', string='相位编码时长 (ms)', val=2, field='SEQ')
        self.addParameter(key='readoutTime', string='读出时长 (ms)', val=3.0, field='SEQ')

        self.addParameter(key='shimming', string='线性匀场 [x, y, z]', val=[210.0, 210.0, 895.0], field='OTH')
        self.addParameter(key='gFactor', string='梯度效率 [x, y, z]', val=[1.0, 1.0, 1.0], field='OTH')
        self.addParameter(key='readPadding', string='读出边距 (μs)', val=1.0, field='OTH')

    def sequenceInfo(self):
        print("============ Averaging Acquisition with Relaxation Enhancement ============")
        print("利用RARE 3D的回波链进行平均的SE序列")

    def sequenceTime(self):  # minutes
        self.sequenceAtributes()
        return self.repetitionTime * self.nPoints[1] * self.nPoints[2] / 6e4

    def sequenceAtributes(self):
        super().sequenceAtributes()  # 把mapVals里的键值对读进self里

        read_points = self.etl*np.product(self.nPoints)*hw.oversamplingFactor
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
        self.axes = {
            'rd': self.axesOrientation[0],
            'ph': self.axesOrientation[1],
            'sl': self.axesOrientation[2]
        }

        print('TE_eff: ', round((self.t_e * self.etl + self.tau) / 1e3), 'ms')
        acqTime = self.readoutTime - 2*self.readPadding  # 读出边距
        self.samplingPeriod = acqTime / self.over_samples

        if not np.product(self.voxel) * self.readoutTime > 0:  # 防止输入过程中出现0
            return 0
        self.readoutAmp = 2e6 * self.gFactor[self.axes['rd']] / self.voxel[0] / hw.gammaB / self.readoutTime 
        phAmp = 1e6 * self.gFactor[self.axes['ph']] / self.voxel[1] / hw.gammaB / (self.phaseTime + self.riseTime)
        slAmp = 1e6 * self.gFactor[self.axes['sl']] / self.voxel[2] / hw.gammaB / (self.phaseTime + self.riseTime)
        print('梯度幅值: ', self.readoutAmp, phAmp, slAmp)
        self.phaseGrads = np.linspace(-phAmp, phAmp, self.nPoints[1])  # 相位编码梯度
        self.sliceGrads = np.linspace(-slAmp, slAmp, self.nPoints[2])  # 选层编码梯度
        if self.nPoints[2] < 2:
            self.sliceGrads = np.array([0])

        self.rfReTime = self.rfExTime
        echoSpacingMin = self.rfReTime + hw.blkTime + 2*self.phaseTime + 6*self.riseTime + self.readoutTime
        print('最小回波间隔: ', echoSpacingMin, 'μs')
        if self.t_e < echoSpacingMin:
            self.t_e = echoSpacingMin + 2
            self.mapVals['echoSpacing'] = self.t_e / 1e3

        # 计算ReadOut predephase梯度大小
        self.ROpreAmp = self.readoutAmp*(self.readoutTime + self.riseTime)/(self.phaseTime + self.riseTime)/2
        if np.abs(self.ROpreAmp) > 1:
            print('ReadOut predephase梯度过大！')
            return 0

    def sequenceRun(self, plotSeq=0, demo=False):
        self.sequenceAtributes()
        self.expt = Experiment(lo_freq=self.mapVals['larmorFreq'], rx_t=self.samplingPeriod)
        self.mapVals['samplingRate'] = self.expt.get_rx_ts()[0]  # 采样间隔
        acq_time = self.mapVals['samplingRate'] * self.over_samples
        print('采样率：%dkHz' % (1e3/self.mapVals['samplingRate']/hw.oversamplingFactor))
        
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
            for ch in range(hw.rx_channels):  # 噪声
                self.rxGate(t_start, acq_time, ch)
            t_start += acq_time + 10

            self.rfRecPulse(t_start, self.rfExTime, self.rfExAmp)  # 激发

            # 读出方向predephase
            t_start += hw.blkTime + self.rfExTime
            gradient(t_start, self.phaseTime, self.ROpreAmp, self.axes['rd']) 
        
            t_start += self.tau - self.rfExTime/2
            for i in range(self.etl):
                t_refocus = t_start + self.t_e*i 
                t_echo = t_refocus + self.tau

                # Refocusing pulse
                self.rfRecPulse(t_refocus - hw.blkTime - self.rfReTime/2, 
                                self.rfReTime, self.rfReAmp, np.pi/2)
                
                # Readout
                t_read = t_echo - self.readoutTime/2
                gradient(t_read - self.riseTime, self.readoutTime, self.readoutAmp, self.axes['rd'])
                for ch in range(hw.rx_channels):
                    self.rxGate(t_read + self.readPadding, acq_time, ch)

                # Slice and Phase encoding
                t_phase = t_read - 3*self.riseTime - self.phaseTime
                gradient(t_phase, self.phaseTime, phase_amp, self.axes['ph'])
                gradient(t_phase, self.phaseTime, slice_amp, self.axes['sl'])

                # phase and slice rewind, no spoil
                t_rewind = t_read + self.readoutTime + self.riseTime + 1
                gradient(t_rewind, self.phaseTime, -phase_amp, self.axes['ph'])
                gradient(t_rewind, self.phaseTime, -slice_amp, self.axes['sl'])

        raw_data = np.zeros((
            self.nPoints[2], self.nPoints[1], hw.rx_channels, (self.etl + 1)*self.over_samples
        ), complex)
        try:
            for i in range(self.nPoints[2]):
                for j in range(self.nPoints[1]):
                    self.iniSequence(20, self.shimming)
                    shot(100, self.sliceGrads[i], self.phaseGrads[j])
                    self.endSequence(100 + self.t_r)
                    self.expt.__del__()
                    self.expt = Experiment(self.mapVals['larmorFreq'], self.mapVals['samplingRate'])
                    if not self.floDict2Exp():
                        print('seq CE %d,%d' % (i, j))
                        return 0
                    rxd, _ = self.expt.run()
                    raw_data[i, j] = list(rxd.values())
                    
            self.mapVals['rawData'] = raw_data
            self.mapVals['samplesPerRead'] = self.over_samples
            return True
        except Exception as e:
            print(e)
            return 0
        finally:
            self.expt.__del__()

    def sequenceAnalysis(self):
        raw_data = self.mapVals['rawData']
        for ch in range(hw.rx_channels):
            data_over = raw_data[:, :, ch, :].reshape(-1, self.over_samples)
            data_deci = np.apply_along_axis(decimate, 1, data_over, hw.oversamplingFactor, ftype='fir', zero_phase=True)
            data_full = data_deci.reshape(self.nPoints[2], self.nPoints[1], -1, self.nPoints[0])
            self.mapVals['mse3d_ch%d' % ch] = data_full[:, :, 1:, :]
            self.mapVals['noise3d_ch%d' % ch] = data_full[:, :, 0, :].reshape(self.nPoints[2], self.nPoints[1], -1)
            ksp = self.mapVals['ksp3d_ch%d' % ch] = np.mean(self.mapVals['mse3d_ch%d' % ch], 2)  # 对回波链取平均
            self.mapVals['img3d_ch%d' % ch] = ifftshift(ifftn(ifftshift(ksp)))  # 重建

        img = self.mapVals['img3d_ch0']    
        abs_img = np.abs(img)
        self.output = self.out = [{
            'widget': 'image',
            'data': np.concatenate((abs_img/np.max(abs_img)*2 - 1, np.angle(img) / np.pi)),
            'xLabel': '相位编码',
            'yLabel': '频率编码',
            'title': '幅值图与相位图',
            'row': 0,
            'col': 0
        }, {
            'widget': 'image',
            'data': np.log(np.abs(self.mapVals['ksp3d_ch0'])),
            'xLabel': 'mV',
            'yLabel': 'ms',
            'title': 'k-Space',
            'row': 0,
            'col': 1
        }]
        self.saveRawData()
        return self.output

    @property
    def over_samples(self):
        return self.nPoints[0] * hw.oversamplingFactor
