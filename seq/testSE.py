"""
@author: T. Guallart Naval
MRILAB @ I3M
"""

import experiment as ex
import numpy as np
import seq.mriBlankSeq as blankSeq  # Import the mriBlankSequence for any new sequence.
import scipy.signal as sig
import configs.hw_config as hw
from plotview.spectrumplot import SpectrumPlot, Spectrum3DPlot
from scipy.optimize import curve_fit
import pyqtgraph as pg
import csv



class testSE(blankSeq.MRIBLANKSEQ):
    def __init__(self):
        super(testSE, self).__init__()
        # Input the parameters
        self.addParameter(key='seqName', string='testSE', val='testSE')
        self.addParameter(key='larmorFreq', string='Larmor frequency (MHz)', val=3.075, field='RF')
        self.addParameter(key='rfExAmp', string='RF excitation amplitude (a.u.)', val=0.3, field='RF')
        self.addParameter(key='rfReAmp', string='RF refocusing amplitude (a.u.)', val=0.3, field='RF')
        self.addParameter(key='rfExTime', string='RF excitation time (us)', val=36.0, field='RF')
        self.addParameter(key='rfReTime', string='RF refocusing time (us)', val=72.0, field='RF')
        self.addParameter(key='phaseRe', string='RF refocusing phase', val=np.pi/2, field='RF')
        self.addParameter(key='echoSpacing', string='Echo spacing (ms)', val=10.0, field='SEQ')
        self.addParameter(key='etl', string='Echo train length', val=1, field='SEQ')
        self.addParameter(key='repetitionTime', string='Repetition time (ms)', val=50., field='SEQ')
        self.addParameter(key='nRepetitions', string='Number of repetitions ', val=60, field='SEQ')
        self.addParameter(key='nScans', string='Number of scans ', val=60, field='SEQ')
        self.addParameter(key='acqCenter', string='Acq center (ms)', val=0.0, field='SEQ')
        self.addParameter(key='nPoints', string='nPoints', val=90, field='IM')
        self.addParameter(key='acqTime', string='Acquisition time (ms)', val=4.0, field='SEQ')
        self.addParameter(key='ttlExtra', string='TTL (1-pi/2 pulse; 2-pi pulse) (ms)', val=2, field='SEQ')
        self.addParameter(key='plotOption', string='Plot (0 = angle; 1 = max Echo)', val=0, field='OTH')

    def sequenceInfo(self):
        print(" ")
        print("Testing SE")
        print("Author: T. Guallart Naval")
        print("mriLab @ i3M, CSIC, Spain")
        print("This sequence runs a spin echo without gradients")

    def sequenceTime(self):
        repetitionTime = self.mapVals['repetitionTime']*1e-3
        nRepetitions  =  self.mapVals['nRepetitions']
        nScans = self.mapVals['nScans']
        return(repetitionTime*nRepetitions*nScans/60)  # minutes, scanTime

    def sequenceRun(self, plotSeq=0):
        init_gpa = False  # Starts the gpa

        seqName = self.mapVals['seqName']
        larmorFreq = self.mapVals['larmorFreq']
        rfExAmp = self.mapVals['rfExAmp']
        rfExTime = self.mapVals['rfExTime']
        rfReAmp = self.mapVals['rfReAmp']
        rfReTime = self.mapVals['rfReTime']
        phaseRe = self.mapVals['phaseRe']
        echoSpacing = self.mapVals['echoSpacing']
        etl = self.mapVals['etl']
        repetitionTime = self.mapVals['repetitionTime']
        nRepetitions  =  self.mapVals['nRepetitions']
        nScans = self.mapVals['nScans']
        nPoints = self.mapVals['nPoints']
        acqTime = self.mapVals['acqTime']
        acqCenter = self.mapVals['acqCenter']
        ttlExtra = self.mapVals['ttlExtra']

        def createSequence():
            # Initialize time
            t0 = 25
            tEx = 20e3

            for nRep in range(nRepetitions):
                # Excitation pulse

                t0 = tEx - hw.blkTime - rfExTime / 2
                t0Ex = t0
                self.rfRecPulse(t0, rfExTime, rfExAmp, 0)
                # self.ttl(t0, rfExTime+hw.blkTime, channel=0)
                # if ttlExtra == 1:
                #     self.ttl(t0, rfExTime, channel=1)
                for nEcho in range(etl):
                    # Refocusing pulse
                    t0 = tEx + nEcho*echoSpacing + echoSpacing/2 - hw.blkTime - rfReTime / 2
                    self.rfRecPulse(t0, rfReTime, rfReAmp, phaseRe*np.pi/180)
                    # self.ttl(t0, rfReTime+hw.blkTime, channel=0)
                    # self.ttl(t0, rfReTime+hw.blkTime, channel=1)
                    # if ttlExtra == 2:
                    #     self.ttl(t0, rfReTime, channel=1)

                    # Rx gate
                    tEcho = tEx + nEcho*echoSpacing + echoSpacing - acqCenter
                    t0 = tEcho - acqTime / 2
                    self.rxGate(t0, acqTime)
                    # self.ttl(t0, acqTime, channel=0)


                # Update time for next repetition
                tEx = tEx + repetitionTime
            self.endSequence(20e3 + repetitionTime*nRepetitions)

        # Time variables in us
        echoSpacing = echoSpacing * 1e3
        repetitionTime = repetitionTime * 1e3
        acqTime = acqTime * 1e3
        acqCenter = acqCenter * 1e3

        # Initialize the experiment
        # bw = 50
        bw = nPoints / acqTime #* hw.oversamplingFactor  # MHz
        samplingPeriod = 1 / bw  # us
        # gpa_fhdo_offset_time= (1 / 0.2 / 3.1)
        self.expt = ex.Experiment(lo_freq=larmorFreq, rx_t=samplingPeriod, init_gpa=init_gpa, gpa_fhdo_offset_time= 0)
        samplingPeriod = self.expt.get_rx_ts()[0]
        bw = 1 / samplingPeriod #/ hw.oversamplingFactor  # MHz
        acqTime = nPoints / bw  # us
        self.mapVals['bw'] = bw
        createSequence()

        if plotSeq == 0:
            # Run the experiment and get data
            print('Running...')
            dataFull = []
            spectrumFull = []
            for nScan in range(nScans):
                rxd, msgs = self.expt.run()
                print(msgs)
                self.mapVals['dataFull'] = rxd['rx0'] #* 13.788
                data = rxd['rx0'] #* 13.788
                # data = sig.decimate(data, hw.oversamplingFactor, ftype='fir', zero_phase=True)
                dataFull = np.concatenate((dataFull, data), axis=0)
                # data = np.reshape(data, (nRepetitions, nPoints))
                # for nRep in range(nRepetitions):
                #     spectrum = np.abs(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(data[nRep]))))
                #     spectrumFull = np.concatenate((spectrumFull, spectrum), axis=0)

            data = np.reshape(dataFull, (nRepetitions*nScans*etl, -1))
            self.mapVals['data'] = data
            spectrum = np.reshape(spectrumFull, (nRepetitions * nScans * etl, -1))
            self.mapVals['spectrum'] = spectrum

        self.expt.__del__()

        return 0

    def sequenceAnalysis(self, obj=''):
        self.saveRawData()
        plotOption = self.mapVals['plotOption']
        data = self.mapVals['data']
        # spectrum = self.mapVals['spectrum']
        # bw = self.mapVals['bw']

        # magnitude = Spectrum3DPlot(np.abs(data), title="Magnitude")
        # magnitudeWidget = magnitude.getImageWidget()
        #
        # phase = Spectrum3DPlot(np.angle(data), title="Phase")
        # phaseWidget = phase.getImageWidget()
        #
        # win = pg.LayoutWidget()
        # win.resize(300, 1000)
        # win.addWidget(magnitudeWidget, row=0, col=0)
        # win.addWidget(phaseWidget, row=0, col=1)
        # return([win])

        # data = np.reshape(data, -1)
        acqTime = self.mapVals['acqTime']
        nRepetitions = self.mapVals['nRepetitions']
        nScans = self.mapVals['nScans']
        nPoints = self.mapVals['nPoints']
        etl = self.mapVals['etl']
        timeVector = np.linspace(0, acqTime*nRepetitions*nScans*etl, num=nPoints*nRepetitions*nScans*etl)
        timeVector = np.transpose(timeVector)


        # # fVector = np.linspace(0, bw*nRepetitions*nScans, nPoints*nRepetitions*nScans)
        # # fVector = np.transpose(fVector)
        # fVectorFull = []
        # fVector = np.linspace(-bw/2, bw/2, nPoints)
        # for nIndex in range(nRepetitions*nScans):
        #     fVectorFull = np.concatenate((fVectorFull, fVector), axis=0)
        # fVector = np.transpose(fVectorFull)
        #
        data = np.reshape(data, -1)
        # spectrum = np.reshape(spectrum, -1)

        # Plot signal versus time
        magPlotWidget = SpectrumPlot(xData=timeVector,
                                yData=[np.abs(data), np.real(data), np.imag(data)],
                                legend=['abs', 'real', 'imag'],
                                xLabel='Time (ms)',
                                yLabel='Signal amplitude (mV)',
                                title='Magnitude')

        # specPlotWidget = SpectrumPlot(xData=fVector,
        #                              yData=[spectrum],
        #                              legend=['abs'],
        #                              xLabel='f (kHz)',
        #                              yLabel='spectrum amplitude (a. u)',
        #                              title='FFT')
        anglePlotWidget = SpectrumPlot(xData=timeVector,
                                      yData=[np.angle(data)],
                                      legend=['abs', 'real', 'imag'],
                                      xLabel='Time (ms)',
                                      yLabel='Phase (rad)',
                                      title='Phase')

        repetitions = np.linspace(1, nRepetitions*nScans*etl, nRepetitions*nScans*etl)
        data = np.reshape(data, (nRepetitions*nScans*etl, -1))
        phase = np.angle(data[:, int(nPoints/2)])
        phasePlotWidget = SpectrumPlot(xData=repetitions,
                                       yData=[np.unwrap(phase)],
                                       legend=[''],
                                       xLabel='Repetition',
                                       yLabel='Phase (rad)',
                                       title='Phase')


        repetitions = np.linspace(1, nRepetitions*nScans*etl, nRepetitions*nScans*etl)
        data = np.reshape(data, (nRepetitions*nScans*etl, -1))
        maxAbs = np.abs(data[:, int(nPoints/2)])
        maxAbsPlotWidget = SpectrumPlot(xData=repetitions,
                                       yData=[maxAbs],
                                       legend=[''],
                                       xLabel='Repetition',
                                       yLabel='Central point (abs)',
                                       title='CPMG')

        if plotOption == 1:
            self.out = [magPlotWidget, maxAbsPlotWidget]
        else:
            self.out = [magPlotWidget, phasePlotWidget]

        return(self.out)