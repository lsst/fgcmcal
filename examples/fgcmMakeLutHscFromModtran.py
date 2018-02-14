"""
HSC-specific overrides for FgcmMakeLut
"""

# Short-code Filter names
config.filterNames = ('g', 'r', 'i', 'i2', 'z', 'y')
# Each filter maps onto a "standard filter".
# Here, both i and i2 map onto i2 which will be the standard
config.stdFilterNames = ('g', 'r', 'i2', 'i2', 'z', 'y')

# Telescope elevation in meters
config.elevation = 4139.0

# Pressure range (in millibar)
config.pmbRange = [820.0, 835.0]
# Number of PMB steps for LUT
config.pmbSteps = 5

# Water vapor range (in mm)
config.pwvRange = [0.1, 3.6]
# Number of PWV steps for LUT
config.pwvSteps = 11

# Ozone range (Dobson)
config.o3Range = [220.0, 310.0]
# Number of O3 steps for LUT
config.o3Steps = 3

# Aerosol Optical Depth (AOD tau) range (unitless)
config.tauRange = [0.002, 0.35]
# Number of tau steps for LUT
config.tauSteps = 11
# Aerosol normalization wavelength (A)
config.lambdaNorm = 7750.0

# Aerosol Optical depth index (alpha) range (unitless)
config.alphaRange = [0.0, 2.0]
# Number of alpha steps for LUT
config.alphaSteps = 9

# Zenith angle range (degrees)
config.zenithRange = [0.0, 70.0]
# Number of zenith angle steps
config.zenithSteps = 21

# Standard value for PMB
config.pmbStd = 828.0
# Standard value for PWV
config.pwvStd = 1.5
# Standard value for O3
config.o3Std = 263.0
# Standard value for tau
config.tauStd = 0.030
# Standard value for alpha
config.alphaStd = 1.0
# Standard value for airmass
config.airmassStd = 1.1
