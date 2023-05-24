import numpy as np

import configs.hw_config as hw
import controller.experiment_gui as ex
import seq.mriBlankSeq as blankSeq


class GRE2Denoise(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(GRE2Denoise, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='梯度回波成像', val='GRE2Denoise')

        self.addParameter(key='nPoints', string='像素点数', val=256, field='IM')
        self.addParameter(key='nScans', string='平均次数', val=1, field='IM')
        self.addParameter(key='phaseAmp', string='相位编码步进 (o.u.)', val=0.001, field='IM')
        self.addParameter(key='readAmp', string='读出梯度幅值 (o.u.)', val=0.1, field='IM')

        self.addParameter(key='flipAngle', string='翻转角（°）', val=30, field='RF')
        self.addParameter(key='larmorFreq', string='中心频率 (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(key='bwCalib', string='调频带宽 (MHz)', val=0.02, field='RF')

        self.addParameter(key='echoSpacing', string='TE (ms)', val=1.2, field='SEQ')
        self.addParameter(key='repetitionTime', string='TR (ms)', val=100, field='SEQ')
        self.addParameter(key='spoilDelay', string='扰相延迟 (μs)', val=100, field='SEQ')
        self.addParameter(key='readPadding', string='读出边距 (μs)', val=10, field='SEQ')

        self.addParameter(key='shimming', string='线性匀场 [x,y,z]', val=[350.0, 370.0, 1030.0], field='OTH')
        self.addParameter(key='axes', string='[读出，相位，选层]', val=[0, 1, 2], field='OTH')
        self.addParameter(key='raiseTime', string='梯度上升时间 (μs)', val=300, field='OTH')
        self.addParameter(key='raiseSteps', string='梯度上升步数', val=20, field='OTH')

    def sequenceInfo(self):
        print("============ 梯度回波成像 ============")

    def sequenceTime(self):
        flip_angle = self.mapVals['flipAngle'] / 90 * np.pi  # 角度转弧度
        rf_amp = round(flip_angle / (100 * hw.b1Efficiency), 6)
        if self.mapVals.get('rfExAmp') is None:
            self.mapVals['rfExAmp'] = 0.05
        if rf_amp != self.mapVals['rfExAmp']:
            self.mapVals['rfExAmp'] = rf_amp
            print("激发功率：", self.mapVals['rfExAmp'], " (a.u.)")
        return self.mapVals['repetitionTime'] * self.mapVals['nPoints'] * self.mapVals['nScans'] / 6e4

    def sequenceRun(self, plot_seq=0):
        print('扫描用时：', self.sequenceTime())

        num_points = self.mapVals['nPoints']
        num_scans = self.mapVals['nScans']
        phase_amp = self.mapVals['phaseAmp']
        read_amp = self.mapVals['readAmp']
        t_e = self.mapVals['echoSpacing'] * 1e3
        t_r = self.mapVals['repetitionTime'] * 1e3
        spoil_delay = self.mapVals['spoilDelay']
        read_padding = self.mapVals['readPadding']

        shimming = np.array(self.mapVals['shimming']) * 1e-4
        axes = self.mapVals['axes']
        rise_time = self.mapVals['raiseTime']
        g_steps = self.mapVals['raiseSteps']
        rf_amp = self.mapVals['rfExAmp']
        bw_calib = self.mapVals['bwCalib']

        rf_time = self.mapVals['rfExTime'] = 100
        if bw_calib > 0:
            hw.larmorFreq = self.freqCalibration(bw_calib, 0.001)
            print("频率校准：", hw.larmorFreq, " (MHz)")
            hw.larmorFreq = self.freqCalibration(bw_calib)
            self.mapVals['larmorFreq'] = hw.larmorFreq

        dephase_refocus_time = t_e - hw.deadTime - rf_time/2 - 3*rise_time
        refocus_time = dephase_refocus_time + rise_time
        dephase_time = refocus_time/2 - rise_time
        print('dephase: ', dephase_time, 'refocus: ', refocus_time)

        acq_time = refocus_time - 2*read_padding + 6000
        sampling_period = acq_time / num_points / hw.oversamplingFactor
        self.expt = ex.Experiment(lo_freq=self.mapVals['larmorFreq'], rx_t=sampling_period)
        self.mapVals['samplingRate'] = self.expt.getSamplingRate()
        acq_time = self.mapVals['samplingRate'] * num_points
        print('采样率：', 1e3/self.mapVals['samplingRate'], ' (kHz)')
        print('采样时长：', acq_time, ' (μs)')

        phase_amp_max = phase_amp * (num_points - 1) / 2
        if phase_amp_max > 1:
            print('相位编码过大！')
            return 0

        def gradient(t, flat, amp, channel):
            self.gradTrap(
                tStart=t,
                gRiseTime=rise_time,
                gFlattopTime=flat,
                gAmp=amp,
                gSteps=g_steps,
                gAxis=channel,
                shimming=shimming
            )

        # --------------------- ↓序列开始↓ ---------------------
        tim = 20
        self.iniSequence(tim, shimming)

        for i in range(num_scans):
            for j in range(num_points):
                tim = 1e5 + (i * num_points + j) * t_r
                self.rfRecPulse(tim, rf_time, rf_amp)

                tim += hw.blkTime + rf_time + hw.deadTime
                gradient(tim, dephase_time, -read_amp, axes[0])
                gradient(tim, dephase_time, (num_points/2 - j) * phase_amp, axes[1])

                tim += dephase_time + 2 * rise_time + 1
                gradient(tim, refocus_time, read_amp, axes[0])

                tim += rise_time + read_padding
                self.rxGateSync(tim, acq_time)
                self.rxGateSync(tim, acq_time, 1)

                tim += refocus_time - read_padding + rise_time + spoil_delay
                for k in range(3):
                    gradient(tim, dephase_time, phase_amp_max, axes[k])

                tim += dephase_time + 2*rise_time + read_padding
                self.rxGateSync(tim, acq_time)
                self.rxGateSync(tim, acq_time, 1)

        self.endSequence(num_points * num_scans * t_r + 2e6)
        # --------------------- ↑序列结束↑ ---------------------

        if not self.floDict2Exp():
            return 0

        if not plot_seq:
            rxd, msg = self.expt.run()
            print(msg)
            self.mapVals['dataOver0'] = rxd['rx0'] * hw.adcFactor
            self.mapVals['dataOver1'] = rxd['rx1'] * hw.adcFactor

        self.expt.__del__()

    def sequenceAnalysis(self):
        num_points = self.mapVals['nPoints']
        num_scans = self.mapVals['nScans']
        data_matrix = [self.mapVals['dataOver0'], self.mapVals['dataOver1']]

        for ch in range(2):
            data_over = np.reshape(data_matrix[ch], (num_scans, -1))
            data_former = np.zeros((num_points, num_points, num_scans)) * 0j
            data_latter = np.zeros((num_points, num_points, num_scans)) * 0j
            for i in range(num_scans):
                data_scan = np.reshape(data_over[i], (num_points, -1))
                n = int(np.size(data_scan[0]) / 2)
                data_over_former = data_scan[:, :n]
                data_over_latter = data_scan[:, n:]
                data_full_former = self.decimate(np.array(data_over_former), num_points)
                data_full_latter = self.decimate(np.array(data_over_latter), num_points)
                data_former[:, :, i] = data_full_former.reshape((num_points, num_points))
                data_latter[:, :, i] = data_full_latter.reshape((num_points, num_points))
            self.mapVals['dataFormer%i' % ch] = data_former
            self.mapVals['dataLatter%i' % ch] = data_latter

        ksp = np.average(self.mapVals['dataFormer0'], 2)
        self.mapVals['img'] = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(ksp)))  # 重建
        self.saveRawData()

        img = np.reshape(self.mapVals['img'], (1, num_points, num_points))
        ksp = np.reshape(ksp, (1, num_points, num_points))
        return [{
            'widget': 'image',
            'data': np.log(np.abs(img)),
            'xLabel': '相位编码',
            'yLabel': '频率编码',
            'title': '幅值图',
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
