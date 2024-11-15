from abc import ABCMeta, abstractmethod
from typing import Literal
import numpy as np
from configs import hw_config as hw, units
from seq import mriBlankSeq as SeqBase
from experiment import Experiment

type LogicalChannel = Literal['slice', 'x', 'y', 'z2']
type NonlinearGradients = dict[
    LogicalChannel, tuple[list[float], list[float]]
]

class NonlinearBase(SeqBase.MRIBLANKSEQ, metaclass=ABCMeta):
    def __init__(self):
        super().__init__()
        self.gradients: NonlinearGradients = {
            'slice': [[], []],
            'x': [[], []],
            'y': [[], []],
            'z2': [[], []]
        }

    def sequenceRun(self, plotSeq=0, demo=False):
        self.sequenceAtributes()
        self.shot()
        
    def sequenceAnalysis(self):
        self.output = []
        self.saveRawData()
        return self.output
    
    @abstractmethod
    def shot(self):
        raise NotImplementedError
    
    def nonlinealTrap(self, t0: float, flatTime: float, slopeTime: float, steps: int, channel: LogicalChannel):
        '''非线性梯度脉冲'''
        self.gradients[channel]

    def logic2physical(self):
        '''4路逻辑梯度转换为40路物理梯度'''
        transformMatrix = np.eye(4, 40)
        np.prod(transformMatrix, self.gradients)
        self.gradTrap(self.gradients['slice'])

class SelectiveExcitationMixin(NonlinearBase, metaclass=ABCMeta):
    def __init__(self):
        super().__init__()
        self.addParameter('larmorFreq', 'Larmor frequency (MHz)', hw.larmorFreq, units.MHz, 'RF')
        self.addParameter('rfExAmp', 'RF excitation amplitude (a.u.)', 0.03, True, 'RF')
        self.addParameter('rfExTime', 'RF excitation time (us)', 200.0, units.us, 'RF')
        self.addParameter('sliceThickness', 'Slice thickness (cm)', 2.0, units.cm, 'IM')

    def addSelectiveExcitation(self, t0: float):
        self.rfSincPulse(t0, self.mapVals['rfExTime'], self.mapVals['rfExAmp'], 0, 3)
        # self.gradTrap(t0 - hw.grad_rise_time, hw.grad_rise_time, self.mapVals['rfExTime'], self.mapVals['sliceThickness'], hw.grad_steps)
        self.nonlinealTrap(t0, self.mapVals['rfExTime'], hw.grad_rise_time, hw.grad_steps, 'slice')
