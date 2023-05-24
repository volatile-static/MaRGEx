import numpy as np
larmorFreq = 4.4286 # MHz
b1Efficiency = np.pi / (0.117 * 100) # rads / (a.u. * us)
gFactor = [0.025, 0.040, 0.035] # (X, Y, Z) in T/m/o.u.
slewRate = 1000 # us/o.u., slew rate for gradient rises
blkTime = 35 # us, blanking time of Barthel's RFPA
gradDelay = 9 # Gradient amplifier delay (us)
deadTime = 160 # us, RF coil dead time
adcFactor = 16 # mV/adcUnit
antenna_dict = {
    "头": np.pi / (0.117 * 100), 
    "膝": np.pi / (0.117 * 100), 
    "腕": np.pi / (0.055 * 100)
}
