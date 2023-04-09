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
            dataFull = self.decimate(rxd['rx0'], nScans)

            matrix = np.reshape(dataFull, (nScans, -1))
            peakVals = []
            for ii in range(nScans):
                freq = np.fft.fftshift(np.fft.fft(
                    np.fft.fftshift(matrix[ii]), n=matrix[ii].size))
                peakVals.append(np.max(np.abs(freq)))
            snr = 20*np.log10(np.mean(peakVals)/np.std(peakVals))
            print('SNR = ', snr, 'dB')
            self.mapVals['snr'] = snr

            # Average data
            data = np.average(matrix, axis=0)
            self.mapVals['data'] = data

            # Save data to sweep plot (single point)
            self.mapVals['sampledPoint'] = data[0]

        self.expt.__del__()

    def sequenceAnalysis(self):
        # Signal and spectrum from 'fir' and decimation
        signal = self.mapVals['data']
        bw = self.mapVals['bw']*1e3  # kHz
        nPoints = self.mapVals['nPoints']
        deadTime = hw.deadTime*1e-3  # ms
        rfExTime = self.mapVals['rfExTime']*1e-3  # ms
        tVector = np.linspace(rfExTime/2 + deadTime + 0.5/bw,
                              rfExTime/2 + deadTime + (nPoints-0.5)/bw, nPoints)
        fVector = np.linspace(-bw/2, bw/2, nPoints)
        spectrum = np.abs(np.fft.ifftshift(
            np.fft.ifftn(np.fft.ifftshift(signal))))
        fitedLarmor = self.mapVals['larmorFreq'] + \
            fVector[np.argmax(np.abs(spectrum))] * 1e-3

        # Get the central frequency
        print('Larmor frequency: %1.5f MHz' % fitedLarmor)
        self.mapVals['signalVStime'] = [tVector, signal]
        self.mapVals['spectrum'] = [fVector, spectrum]

        self.saveRawData()

        # Add time signal to the layout
        result1 = {'widget': 'curve',
                   'xData': tVector,
                   'yData': [np.abs(signal), np.real(signal), np.imag(signal)],
                   'xLabel': 'Time (ms)',
                   'yLabel': 'Signal amplitude (mV)',
                   'title': 'Signal vs time',
                   'legend': ['abs', 'real', 'imag'],
                   'row': 0,
                   'col': 0}

        # Add frequency spectrum to the layout
        result2 = {'widget': 'curve',
                   'xData': fVector,
                   'yData': [spectrum],
                   'xLabel': 'Frequency (kHz)',
                   'yLabel': 'Spectrum amplitude (a.u.)',
                   'title': 'Spectrum',
                   'legend': [''],
                   'row': 1,
                   'col': 0}

        # create self.out to run in iterative mode
        self.out = [result1, result2]

        return self.out
