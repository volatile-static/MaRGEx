import numpy as np

import configs.hw_config as hw
import controller.experiment_gui as ex
import seq.mriBlankSeq as blankSeq


class GRE2D5(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(GRE2D5, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='梯度回波成像', val='GRE2D')

        self.addParameter(key='nPoints', string='像素点数', val=256, field='IM')
        self.addParameter(key='nScans', string='平均次数', val=1, field='IM')
        self.addParameter(key='phaseAmp', string='相位编码步进 (o.u.)', val=0.00001, field='IM')
        self.addParameter(key='readAmp', string='读出梯度幅值 (o.u.)', val=0.001, field='IM')
        self.addParameter(key='sliceAmp', string='选层梯度幅值 (o.u.)', val=0.001, field='IM')

        self.addParameter(key='rfExAmp', string='激发功率 (a.u.)', val=0.08, field='RF')
        self.addParameter(key='rfExTime', string='激发时长 (μs)', val=500.0, field='RF')
        self.addParameter(key='larmorFreq', string='中心频率 (MHz)', val=hw.larmorFreq, field='RF')
        self.addParameter(key='bwCalib', string='调频带宽 (MHz)', val=0.02, field='RF')

        self.addParameter(key='echoSpacing', string='TE (ms)', val=1.2, field='SEQ')
        self.addParameter(key='repetitionTime', string='TR (ms)', val=100, field='SEQ')
        self.addParameter(key='addReadTime', string='额外读出时长 (ms)', val=6.0, field='SEQ')
        self.addParameter(key='readPadding', string='读出边距 (μs)', val=10, field='SEQ')
        self.addParameter(key='spoilDelay', string='扰相延迟 (μs)', val=100, field='SEQ')

        self.addParameter(key='shimming', string='线性匀场 [x,y,z]', val=[210.0, 210.0, 895.0], field='OTH')
        self.addParameter(key='axes', string='[读出，相位，选层]', val=[0, 1, 2], field='OTH')
        self.addParameter(key='raiseTime', string='梯度上升时间 (μs)', val=300, field='OTH')
        self.addParameter(key='raiseSteps', string='梯度上升步数', val=20, field='OTH')

    def sequenceInfo(self):
        print("============ 梯度回波成像 ============")

    def sequenceTime(self):
        return self.mapVals['repetitionTime'] * self.mapVals['nPoints'] * self.mapVals['nScans'] / 6e4

    def sequenceRun(self, plot_seq=0):
        print('扫描用时：', self.sequenceTime())

        num_points = self.mapVals['nPoints']
        num_scans = self.mapVals['nScans']
        phase_amp = self.mapVals['phaseAmp']
        read_amp = self.mapVals['readAmp']
        slice_amp = self.mapVals['sliceAmp']
        t_e = self.mapVals['echoSpacing'] * 1e3
        t_r = self.mapVals['repetitionTime'] * 1e3
        spoil_delay = self.mapVals['spoilDelay']
        read_padding = self.mapVals['readPadding']
        add_read_time = self.mapVals['addReadTime'] * 1e3

        shimming = np.array(self.mapVals['shimming']) * 1e-4
        axes = self.mapVals['axes']
        rise_time = self.mapVals['raiseTime']
        g_steps = self.mapVals['raiseSteps']
        rf_amp = self.mapVals['rfExAmp']
        rf_time = self.mapVals['rfExTime']
        bw_calib = self.mapVals['bwCalib']

        if bw_calib > 0:
            hw.larmorFreq = self.freqCalibration(bw_calib, 0.001)
            print("频率校准：", hw.larmorFreq, " (MHz)")
            hw.larmorFreq = self.freqCalibration(bw_calib)
            self.mapVals['larmorFreq'] = hw.larmorFreq

        dephase_refocus_time = t_e - hw.deadTime - rf_time/2 - 3*rise_time
        refocus_time = dephase_refocus_time + rise_time
        dephase_time = refocus_time/2 - rise_time
        select_amp = slice_amp*(rf_time + 2*rise_time)/(dephase_time + 2*rise_time)

        acq_time = refocus_time - 2*read_padding + add_read_time
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
                self.rfSincPulse(tim, rf_time, rf_amp, 0, 3)
                gradient(tim + hw.blkTime - rise_time, rf_time, slice_amp, axes[2])

                tim += hw.blkTime + rf_time + rise_time + 1
                gradient(tim, dephase_time, -read_amp, axes[0])
                gradient(tim, dephase_time, (num_points/2 - j) * phase_amp, axes[1])
                gradient(tim, dephase_time, -select_amp, axes[2])

                tim += dephase_time + 2 * rise_time + 1
                gradient(tim, refocus_time, read_amp, axes[0])

                tim += rise_time + read_padding
                self.rxGateSync(tim, acq_time)

                tim += refocus_time - read_padding + rise_time + spoil_delay
                for k in range(3):
                    gradient(tim, dephase_time, phase_amp_max, axes[k])

        self.endSequence(num_points * num_scans * t_r + 2e6)
        # --------------------- ↑序列结束↑ ---------------------

        if not self.floDict2Exp():
            return 0

        if not plot_seq:
            rxd, msg = self.expt.run()
            print(msg)
            self.mapVals['dataOver'] = rxd['rx0'] * hw.adcFactor

        self.expt.__del__()

    def sequenceAnalysis(self):
        num_points = self.mapVals['nPoints']
        data_over = np.reshape(self.mapVals['dataOver'], (self.mapVals['nScans'], -1))
        print('数据维度：', data_over.shape)
        data_average = np.average(data_over, axis=0)
        data_full = self.decimate(data_average, num_points)
        ksp = self.mapVals['ksp'] = np.reshape(data_full, (num_points, num_points))
        self.mapVals['img'] = np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(ksp)))  # 重建
        self.saveRawData()

        img = np.reshape(self.mapVals['img'], (1, num_points, num_points))
        ksp = np.reshape(ksp, (1, num_points, num_points))

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
