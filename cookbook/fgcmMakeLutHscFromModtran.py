"""
HSC-specific overrides for FgcmMakeLut
"""

# Short-code Filter names
config.filterNames = ('g', 'r', 'i', 'z', 'y')
# Each filter maps onto a "standard filter".
# Here, both i and i2 map onto i2 which will be the standard
config.stdFilterNames = ('g', 'r', 'i', 'z', 'y')

# Telescope elevation in meters
config.parameters.elevation = 4139.0

# Pressure range (in millibar)
config.parameters.pmbRange = [820.0, 835.0]
# Number of PMB steps for LUT
config.parameters.pmbSteps = 5

# Water vapor range (in mm)
config.parameters.pwvRange = [0.1, 3.6]
# Number of PWV steps for LUT
config.parameters.pwvSteps = 11

# Ozone range (Dobson)
config.parameters.o3Range = [220.0, 310.0]
# Number of O3 steps for LUT
config.parameters.o3Steps = 3

# Aerosol Optical Depth (AOD tau) range (unitless)
config.parameters.tauRange = [0.002, 0.35]
# Number of tau steps for LUT
config.parameters.tauSteps = 11
# Aerosol normalization wavelength (A)
config.parameters.lambdaNorm = 7750.0

# Aerosol Optical depth index (alpha) range (unitless)
config.parameters.alphaRange = [0.0, 2.0]
# Number of alpha steps for LUT
config.parameters.alphaSteps = 9

# Zenith angle range (degrees)
config.parameters.zenithRange = [0.0, 70.0]
# Number of zenith angle steps
config.parameters.zenithSteps = 21

# Wavelength step size (nm)
config.parameters.lambdaStep = 0.5
# Wavelength range (A)
config.parameters.lambdaRange = [3000.0, 11000.0]

# Standard value for PMB
config.parameters.pmbStd = 828.0
# Standard value for PWV
config.parameters.pwvStd = 1.5
# Standard value for O3
config.parameters.o3Std = 263.0
# Standard value for tau
config.parameters.tauStd = 0.030
# Standard value for alpha
config.parameters.alphaStd = 1.0
# Standard value for airmass
config.parameters.airmassStd = 1.1
