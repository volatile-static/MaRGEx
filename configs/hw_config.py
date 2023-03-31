# Config file for Physio MRI scanner at MRILab, i3M, CSIC, Spain.

# Note: I write Ocra1 Units as o.u.
# Ocra1 gain = 10 V/o.u.
# AE Techrom transductance 5 A/V
# From Ocra1 to current: 50 A/o.u.
# X axis: 25 mT/m/o.u., 0.5 mT/m/A, 2.5 mT/m/V
# Y axis: 40 mT/m/o.u., 0.8 mT/m/A, 4.0 mT/m/V
# Z axis: 35 mT/m/o.u., 0.7 mT/m/A, 3.5 mT/m/V

from configs import tabletop_cfg
# from configs import bedside_cfg

gammaB = 42.56e6 # Hz/T, Gyromagnetic ratio
oversamplingFactor = 6 # Rx oversampling
maxRdPoints = 2**18 # Maximum number of points to be acquired by the red pitaya
maxOrders = 2**14 # Maximum number of orders to be processed by the red pitaya
deadTime = 60 # us, RF coil dead time
cicDelayPoints = 3 # to account for signal delay from red pitaya due to cic filter
addRdPoints = 10 # to account for wrong first points after decimation
adcFactor = 13.788 # mV/adcUnit
arduinoPort = 'COM7'
scanner_name = "Physio V1.01"
fov = [20.0, 20.0, 20.0]
dfov = [0.0, 0.0, 0.0]
bash_path = "D:\Git\git-bash.exe" # use "gnome-terminal" for genome linux
rp_ip_address = "192.168.1.101"
rp_version = "rp-122"
