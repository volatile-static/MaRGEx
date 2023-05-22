"""
Created on Thu June 2 2022
@author: J.M. Algarín, MRILab, i3M, CSIC, Valencia
@email: josalggui@i3m.upv.es
@Summary: All sequences on the GUI must be here
"""

import seq.rare as rare
import seq.rareProtocols as rareProtocols
# import seq.haste as haste
import seq.gre3d as gre
# import seq.petra as petra
import seq.fid as fid
import seq.rabiFlops as rabiFlops
import seq.B1calibration as B1calibration
import seq.cpmg as cpmg
import seq.eddycurrents as eddycurrents
import seq.larmor as larmor
import seq.inversionRecovery as inversionRecovery
import seq.ADCdelayTest as ADCdelayTest
import seq.noise as noise
import seq.shimmingSweep as shimming
import seq.testSE as testSE
# import seq.sliceSelection as sliceSelection
import seq.sweepImage as sweep
# import seq.autoTuning as autoTuning
import seq.localizer as localizer
import seq.softFID as softFID
import seq.gre2denoise as gre2denoise
import seq.gre2d as gre2d
import seq.gre1d as gre1d
import seq.rfopt as rfopt
import seq.snr as snr
import seq.T1SE as t1se

class SNR(snr.SNR):
    def __init__(self): super(SNR, self).__init__()

class SoftFID(softFID.SoftFID):
    def __init__(self): super(SoftFID, self).__init__()

class T1SE(t1se.T1SE):
    def __init__(self): super(T1SE, self).__init__()

class GRE2Denoise(gre2denoise.GRE2Denoise):
    def __init__(self): super(GRE2Denoise, self).__init__()

class GRE2D(gre2d.GRE2D):
    def __init__(self): super(GRE2D, self).__init__()

class GRE1D(gre1d.GRE1D):
    def __init__(self): super(GRE1D, self).__init__()

class RFOPT(rfopt.RFOPT):
    def __init__(self): super(RFOPT, self).__init__()

class RARE(rare.RARE):
    def __init__(self): super(RARE, self).__init__()

class RAREProtocols(rareProtocols.RAREProtocols):
    def __init__(self): super(RAREProtocols, self).__init__()

class testSE(testSE.testSE):
    def __init__(self): super(testSE, self).__init__()

class GRE3D(gre.GRE3D):
    def __init__(self): super(GRE3D, self).__init__()

# class PETRA(petra.PETRA):
#     def __init__(self): super(PETRA, self).__init__()

# class HASTE(haste.HASTE):
#     def __init__(self): super(HASTE, self).__init__()

class FID(fid.FID):
    def __init__(self): super(FID, self).__init__()

class RabiFlops(rabiFlops.RabiFlops):
    def __init__(self): super(RabiFlops, self).__init__()

class B1calibration(B1calibration.B1calibration):
    def __init__(self): super(B1calibration, self).__init__()

class Larmor(larmor.Larmor):
    def __init__(self): super(Larmor, self).__init__()

class Noise(noise.Noise):
    def __init__(self): super(Noise, self).__init__()

class CPMG(cpmg.CPMG):
    def __init__(self): super(CPMG, self).__init__()

class EDDYCURRENTS(eddycurrents.EDDYCURRENTS):
    def __init__(self): super(EDDYCURRENTS, self).__init__()

class IR(inversionRecovery.InversionRecovery):
    def __init__(self): super(IR, self).__init__()

class ADCtest(ADCdelayTest.ADCdelayTest):
    def __init__(self): super(ADCtest, self).__init__()

class Shimming(shimming.ShimmingSweep):
    def __init__(self): super(Shimming, self).__init__()

# class SliceSelection(sliceSelection.SliceSelection):
#     def __init__(self): super(SliceSelection, self).__init__()

class SWEEP(sweep.SweepImage):
    def __init__(self): super(SWEEP, self).__init__()
#
# class AutoTuning(autoTuning.AutoTuning):
#     def __init__(self): super(AutoTuning, self).__init__()

class Localizer(localizer.Localizer):
    def __init__(self): super(Localizer, self).__init__()

"""
Definition of default sequences
"""
defaultsequences = {
    'T1SE': T1SE(),
    'RARE': RARE(),
    'RARE Protocols': RAREProtocols(),
    'GRE2Denoise': GRE2Denoise(),
    'GRE3D': GRE3D(),
    'GRE2D': GRE2D(),
    'GRE1D': GRE1D(),
    # 'PETRA': PETRA(),
    # 'HASTE': HASTE(),
    # 'AutoTuning': AutoTuning(),
    # 'SliceSelection': SliceSelection(),
    'RFOPT': RFOPT(),
    'SoftFID': SoftFID(),
    'FID': FID(),
    'SNR': SNR(),
    'Larmor': Larmor(),
    'Noise': Noise(),
    'RabiFlops': RabiFlops(),
    'B1calibration': B1calibration(),
    'CPMG': CPMG(),
    'EDDYCURRENTS': EDDYCURRENTS(),
    'InversionRecovery': IR(),
    'ADCtest': ADCtest(),
    'Shimming': Shimming(),
    'SWEEP': SWEEP(),
    'testSE': testSE(),
    'Localizer': Localizer(),
}