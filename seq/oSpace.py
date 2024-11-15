import numpy as np
from experiment import Experiment
from configs import units, hw_config as hw
from seq.mriBlankSeq import MRIBLANKSEQ
from seq.nonlinearBase import NonlinearBase, SelectiveExcitationMixin

class OSpace(MRIBLANKSEQ):
    def __init__(self):
        super().__init__()
        self.logicGrad = []
        self.addParameter('seqName', 'OSpaceInfo', 'O-Space')

        self.addParameter('larmorFreq', 'Larmor frequency (MHz)', hw.larmorFreq, units.MHz, 'RF')
        self.addParameter('rfExAmp', 'RF excitation amplitude (a.u.)', .3, True, 'RF')
        self.addParameter('rfExTime', 'RF excitation time (us)', 30., units.us, 'RF')
        self.addParameter('rfLobes', 'Lobes of sinc', 0, True, 'RF')

        self.addParameter('echoTime', 'TE (ms)', .9, units.ms, 'SEQ')
        self.addParameter('repetitionTime', 'TR (ms)', 1000., units.ms, 'SEQ')
        self.addParameter('rxBandwidth', 'Rx bandwidth (kHz)', 20., units.kHz, 'SEQ')
        self.addParameter('r_cp', 'Center point radius (cm)', 12.8, units.cm, 'SEQ')

        self.addParameter('nPoints', 'Number of pixels', [32, 16], True, 'IM')
        self.addParameter('fov', 'Field of view (mm)', [256, 256], units.mm, 'IM')
        self.addParameter('gz2', 'Gz2 (Hz/cm²)', 1., 1 / units.cm**2, 'IM')
        self.addParameter('spoiler', 'Spoiler amplitude (a.u.)', .01, True, 'OTH')

    def sequenceInfo(self):
        print('o-space')

    def sequenceTime(self):
        return self.mapVals['repetitionTime']*units.ms*self.mapVals['nPoints'][1]/60
    
    def sequenceRun(self, plotSeq=0, demo=False):
        def readoutGrad(te, len, amp, ch):
            self.addGrad(te - len - hw.grad_rise_time, len/2 - hw.grad_rise_time, -amp, ch)
            self.addGrad(te - len/2, len, amp, ch)
            # self.gradTrap(te - len - hw.grad_rise_time*2, hw.grad_rise_time, len/2 - hw.grad_rise_time, amp, hw.grad_steps, axis, [0, 0, 0])
            # self.gradTrap(te - len/2 - hw.grad_rise_time, hw.grad_rise_time, len, -amp, hw.grad_steps, axis, [0, 0, 0])

        if not plotSeq:
            self.expt = Experiment(self.larmorFreq, 1e6/self.rxBandwidth)
            samplingRate = 1e6 / self.expt.get_rx_ts()[0]
        else:
            samplingRate = 1 / self.rxBandwidth
        
        readoutTime = self.nPoints[0] / samplingRate
        print(f'{readoutTime=}')

        if self.gz2 > 0:
            gz2 = self.gz2
            gMax = self.gz2 * self.r_cp / hw.gammaB
            fov = samplingRate / gMax / hw.gammaB
            print(f'{fov=}')
        else:
            gMax = samplingRate / self.fov[0] / hw.gammaB  # T/m
            gz2 = gMax / self.r_cp * hw.gammaB  # Hz/m²
        
        self.iniSequence(20, [0, 0, 0])
        t0 = 20
        for θ in range(0, 2*np.pi, self.nPoints[1]):
            gx = np.cos(θ) * self.r_cp * gz2 / hw.gammaB  # T/m
            gy = np.sin(θ) * self.r_cp * gz2 / hw.gammaB  # T/m
            gz = self.fov[1]            
            t_echo = t0 + self.rfExTime/2 + hw.grad_rise_time + self.echoTime
            
            if self.rfLobes > 0:
                self.rfSincPulse(t0 + hw.grad_rise_time, self.rfExAmp, self.rfExTime, 0, self.rfLobes)
                self.addGrad(t0, self.rfExTime, gz, 2)
                self.addGrad(t0 + self.rfExTime + 2*hw.grad_rise_time, self.rfExTime, gz/2, 2)
            else:
                self.rfRecPulse(t0 + hw.grad_rise_time, self.rfExAmp, self.rfExTime)

            self.addGrad(t_echo + readoutTime/2 + hw.grad_rise_time, readoutTime/2, self.spoiler, 2)
            readoutGrad(t_echo, readoutTime, gx, 0)
            readoutGrad(t_echo, readoutTime, gy, 1)
            readoutGrad(t_echo, readoutTime, gz2, 3)

            self.rxGate(t_echo - readoutTime/2, readoutTime, 0)
            t0 += self.repetitionTime

        self.logical2physical()
        self.endSequence(t0)

        if not plotSeq and self.floDict2Exp():
            rxd, _ = self.expt.run(demo)
            self.mapVals['rawData'] = np.reshape(rxd['rx0'] * hw.adcFactor, self.nPoints)
        self.mapVals['readoutTime'] = readoutTime
        self.mapVals['samplingRate'] = samplingRate
        self.mapVals['gz2'] = gz2
        return True

    def sequenceAnalysis(self, mode=None):
        self.output = []
        return self.output
    
    def addGrad(self, t0, len, amp, ch):
        self.logicGrad.append({
            'ch': ch,
            't0': t0,
            'len': len,
            'amp': amp
        })
    
    def logical2physical(self):
        mat = np.eye(4, 40)
        gradSeries = {}
        for d in self.logicGrad:
            k = (d['t0'], d['len'])
            if not k in gradSeries:
                gradSeries[k] = np.zeros(4)
            gradSeries[k][d['ch']] = d['amp']
        for (t0, len), amps in gradSeries.items():
            phy = np.matmul(amps, mat)
            for i in enumerate(phy):
                self.gradTrap(t0, hw.grad_rise_time, len, phy[i], hw.grad_steps, i, [0, 0, 0])

# class OSpace(SelectiveExcitationMixin, NonlinearBase):
#     def __init__(self):
#         super().__init__()
#         self.addParameter('seqName', 'Sequence name', 'OSpace', 'OSpace')

#     def sequenceInfo(self):
#         print("============ O-Space ============")

#     def sequenceTime(self):
#         self.sequenceAtributes()
#         return 0
    
#     def shot():
#         pass
