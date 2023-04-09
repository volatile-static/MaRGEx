import numpy as np
larmorFreq = 20.82 # MHz
b1Efficiency = np.pi / (0.05 * 50) # rads / (a.u. * us)
gFactor = [2.8, 2.8, 3.2] # (X, Y, Z) in T/m/o.u.
slewRate = 1000 # us/o.u., slew rate for gradient rises
stepsRate = 200 # steps/o.u., steps rate for gradient rises
blkTime = 35 # us, blanking time of Barthel's RFPA
gradDelay = 9 # Gradient amplifier delay (us)
antenna_dict = {"收发一体": np.pi / (0.05 * 50), "收发分离": np.pi / (0.05 * 50)}
