import numpy as np
larmorFreq = 63.8 # MHz
b1Efficiency = np.pi / (0.117 * 100) # rads / (a.u. * us)
gFactor = [0.025, 0.040, 0.035] # (X, Y, Z) in T/m/o.u.
slewRate = 1000 # us/o.u., slew rate for gradient rises
blkTime = 35 # us, blanking time of Barthel's RFPA
gradDelay = 9 # Gradient amplifier delay (us)
deadTime = 160 # us, RF coil dead time
adcFactor = 16 # mV/adcUnit
gateActiveLow = True
antenna_dict = {
    "共轭": np.pi / (0.117 * 100), 
}

bash_path = "gnome-terminal"
rp_ip_address = "192.168.5.20"
