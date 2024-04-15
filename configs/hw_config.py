# Config file for Physio MRI scanner at MRILab, i3M, CSIC, Spain.

# Note: I write Ocra1 Units as o.u.
# Ocra1 gain = 10 V/o.u.
# AE Techrom transductance 5 A/V
# From Ocra1 to current: 50 A/o.u.
# X axis: 25 mT/m/o.u., 0.5 mT/m/A, 2.5 mT/m/V
# Y axis: 40 mT/m/o.u., 0.8 mT/m/A, 4.0 mT/m/V
# Z axis: 35 mT/m/o.u., 0.7 mT/m/A, 3.5 mT/m/V
import numpy as np

gFactor = [1.4, 1.4, 1.6] # (X, Y, Z) in T/m/o.u.
grad_rise_time = 400e-6 # s, time for gradient ramps
grad_steps = 16 # steps to gradient ramps
slewRate = 300 # us/o.u.
gammaB = 42.56e6 # Hz/T, Gyromagnetic ratio
blkTime = 35 # us, blanking time of Barthel's RFPA
gradDelay = 9 # Gradient amplifier delay (us)
oversamplingFactor = 6 # Rx oversampling
maxRdPoints = 2**18 # Maximum number of points to be acquired by the red pitaya
maxOrders = 2**14 # Maximum number of orders to be processed by the red pitaya
deadTime = 60 # us, RF coil dead time
b1Efficiency = np.pi / (0.125 * 50) # rads / (a.u. * us)
larmorFreq = 20.798 # MHz
cic_delay_points = 3 # to account for signal delay from red pitaya due to cic filter
addRdPoints = 10 # to account for wrong first points after decimation
adcFactor = 100 # mV/adcUnit
scanner_name = "TableTop 0.5T"
antenna_dict = {"RF01": np.pi / (0.125 * 50), "RF02": np.pi / (0.125 * 50)}
fov = [20.0, 20.0, 20.0]
dfov = [0.0, 0.0, 0.0]
bash_path = "D:\Archivos de Programa\Git\git-bash.exe" # use "gnome-terminal" for genome linux
rp_ip_address = "192.168.5.20"
rp_version = "rp-122"
lnaGain = 50 # dB
temperature = 293 # k
shimming_factor = 1e-5

# Arduinos
ard_sn_autotuning = '242353133363518050E0'
ard_sn_interlock = '242353133363518050E1'
