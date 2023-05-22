import seq.mriBlankSeq as blankSeq
import configs.hw_config as hw
import numpy as np
import controller.experiment_gui as ex
import os
import sys
# *****************************************************************************
# Add path to the working directory
path = os.path.realpath(__file__)
ii = 0
for char in path:
    if (char == '\\' or char == '/') and path[ii+1:ii+14] == 'PhysioMRI_GUI':
        sys.path.append(path[0:ii+1]+'PhysioMRI_GUI')
        sys.path.append(path[0:ii+1]+'marcos_client')
    ii += 1
# ******************************************************************************


class SNR(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(SNR, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='SNRinfo', val='SNR')
        self.addParameter(
            key='larmorFreq', string='Larmor frequency (MHz)', val=20.8, field='RF')
        self.addParameter(
            key='rfExAmp', string='RF excitation amplitude (a.u.)', val=0.03, field='RF')
        self.addParameter(
            key='rfExTime', string='RF excitation time (us)', val=50.0, field='RF')
        self.addParameter(
            key='nScans', string='Number of scans', val=10, field='SEQ')
        self.addParameter(key='repetitionTime',
                          string='Repetition time (ms)', val=100., field='SEQ')
        self.addParameter(
            key='acqTime', string='Acquisition time (ms)', val=1.0, field='SEQ')
        self.addParameter(
            key='nPoints', string='Number of points', val=100, field='IM')
        self.addParameter(key='shimming', string='Shimming (*1e4)',
                          val=[0, 0, 666], field='OTH')

    def sequenceInfo(self):
        print("多次采集FID计算信噪比")

    def sequenceTime(self):
        nScans = self.mapVals['nScans']
        repetitionTime = self.mapVals['repetitionTime']*1e-3
        return (repetitionTime*nScans/60)  # minutes, scanTime

    def sequenceRun(self, plotSeq=0):
        # Create input parameters
        nScans = self.mapVals['nScans']
        larmorFreq = self.mapVals['larmorFreq']  # MHz
        rfExAmp = self.mapVals['rfExAmp']
        rfExTime = self.mapVals['rfExTime']  # us
        repetitionTime = self.mapVals['repetitionTime']*1e3  # us
        acqTime = self.mapVals['acqTime']*1e3  # us
        nPoints = self.mapVals['nPoints']
        shimming = np.array(self.mapVals['shimming'])*1e-4

        # Miscellaneus
        bw = nPoints/acqTime  # MHz

        def createSequence():
            # Shimming
            # shimming is turned on 20 us after experiment beginning
            self.iniSequence(20, shimming)
            self.rfRecPulse(2000, rfExTime, rfExAmp)

            for scan in range(nScans):
                tEx = 2000 + repetitionTime*(scan + 1) + hw.blkTime + rfExTime / 2

                # Excitation pulse
                t0 = tEx - hw.blkTime - rfExTime / 2
                self.rfRecPulse(t0, rfExTime, rfExAmp)

                # Rx gate
                t0 = tEx + rfExTime / 2 + hw.deadTime
                self.rxGateSync(t0, acqTime)
                self.rxGateSync(t0, acqTime, 1)

            self.endSequence(repetitionTime*(nScans + 1))

        # Initialize the experiment
        samplingPeriod = 1 / bw  # us
        self.expt = ex.Experiment(lo_freq=larmorFreq, rx_t=samplingPeriod)
        samplingPeriod = self.expt.getSamplingRate()
        bw = 1 / samplingPeriod
        acqTime = nPoints / bw  # us
        self.mapVals['acqTime'] = acqTime*1e-3  # ms
        self.mapVals['bw'] = bw  # MHz
        createSequence()
        if self.floDict2Exp():
            pass
        else:
            return 0

        if not plotSeq:
            # Run the experiment and get data
            rxd, msgs = self.expt.run()

            # Decimate the signal
            dataFull0 = self.decimate(rxd['rx0'], nScans)
            dataFull1 = self.decimate(rxd['rx1'], nScans)
            matrix1 = np.reshape(dataFull1, (nScans, -1))
            matrix0 = np.reshape(dataFull0, (nScans, -1))
            peakVals = []
            for ii in range(nScans):
                freq = np.fft.fftshift(np.fft.fft(
                    np.fft.fftshift(matrix0[ii]), n=matrix0[ii].size))
                peakVals.append(np.max(np.abs(freq)))
            snr = 20*np.log10(np.mean(peakVals)/np.std(peakVals))
            print('SNR = ', snr, 'dB')
            self.mapVals['snr'] = snr

            # Average data
            data0 = np.average(matrix0, axis=0)
            data1 = np.average(matrix1, axis=0)
            self.mapVals['data0'] = data0
            self.mapVals['data1'] = data1

        self.expt.__del__()

    def sequenceAnalysis(self):
        # Signal and spectrum from 'fir' and decimation
        signal0 = self.mapVals['data0']
        signal1 = self.mapVals['data1']
        bw = self.mapVals['bw']*1e3  # kHz
        nPoints = self.mapVals['nPoints']
        deadTime = hw.deadTime*1e-3  # ms
        rfExTime = self.mapVals['rfExTime']*1e-3  # ms
        tVector = np.linspace(rfExTime/2 + deadTime + 0.5/bw,
                              rfExTime/2 + deadTime + (nPoints-0.5)/bw, nPoints)
        fVector = np.linspace(-bw/2, bw/2, nPoints)
        spectrum0 = np.abs(np.fft.ifftshift(
            np.fft.ifftn(np.fft.ifftshift(signal0))))
        spectrum1 = np.abs(np.fft.ifftshift(
            np.fft.ifftn(np.fft.ifftshift(signal0))))
        fitedLarmor = self.mapVals['larmorFreq'] + \
            fVector[np.argmax(np.abs(spectrum0))] * 1e-3

        # Get the central frequency
        print('Larmor frequency: %1.5f MHz' % fitedLarmor)
        self.mapVals['signalVStime'] = [tVector, signal0]
        self.mapVals['spectrum'] = [fVector, spectrum0]

        self.saveRawData()

        # Add time signal to the layout
        result1 = {'widget': 'curve',
                   'xData': tVector,
                   'yData': [np.abs(signal0), np.abs(signal1)],
                   'xLabel': 'Time (ms)',
                   'yLabel': 'Signal amplitude (mV)',
                   'title': 'Signal vs time',
                   'legend': ['abs0', 'abs1'],
                   'row': 0,
                   'col': 0}

        # Add frequency spectrum to the layout
        result2 = {'widget': 'curve',
                   'xData': fVector,
                   'yData': [spectrum0, spectrum1],
                   'xLabel': 'Frequency (kHz)',
                   'yLabel': 'Spectrum amplitude (a.u.)',
                   'title': 'Spectrum',
                   'legend': ['0', '1'],
                   'row': 1,
                   'col': 0}

        # create self.out to run in iterative mode
        self.out = [result1, result2]

        return self.out
