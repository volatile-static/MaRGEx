# Config file for Physio MRI scanner at MRILab, i3M, CSIC, Spain.

# Note: I write Ocra1 Units as o.u.
# Ocra1 gain = 10 V/o.u.
# AE Techrom transductance 5 A/V
# From Ocra1 to current: 50 A/o.u.
# X axis: 25 mT/m/o.u., 0.5 mT/m/A, 2.5 mT/m/V
# Y axis: 40 mT/m/o.u., 0.8 mT/m/A, 4.0 mT/m/V
# Z axis: 35 mT/m/o.u., 0.7 mT/m/A, 3.5 mT/m/V

gammaB = 42.56e6 # Hz/T, Gyromagnetic ratio
oversamplingFactor = 6 # Rx oversampling
maxRdPoints = 2**18 # Maximum number of points to be acquired by the red pitaya
maxOrders = 2**14 # Maximum number of orders to be processed by the red pitaya
cicDelayPoints = 3 # to account for signal delay from red pitaya due to cic filter
addRdPoints = 10 # to account for wrong first points after decimation
arduinoPort = 'COM7'
scanner_name = "Physio V1.01"
fov = [20.0, 20.0, 20.0]
dfov = [0.0, 0.0, 0.0]
bash_path =  "gnome-terminal"  #"D:\Git\git-bash.exe"
rp_ip_address = "192.168.1.103"#"10.14.102.162"
rp_version = "rp-122"

# from configs import tabletop_cfg as cfg
from configs import bedside_cfg as cfg

larmorFreq = cfg.larmorFreq # MHz
b1Efficiency = cfg.b1Efficiency # rads / (a.u. * us)
gFactor = cfg.gFactor # (X, Y, Z) in T/m/o.u.
slewRate = cfg.slewRate # us/o.u., slew rate for gradient rises
stepsRate = cfg.stepsRate # steps/o.u., steps rate for gradient rises
blkTime = cfg.blkTime # us, blanking time of Barthel's RFPA
gradDelay = cfg.gradDelay # Gradient amplifier delay (us)
antenna_dict = cfg.antenna_dict
deadTime = cfg.deadTime
adcFactor = cfg.adcFactor
