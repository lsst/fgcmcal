"""
HSC-specific overrides for FgcmMakeLut
"""

# Short-code Filter names
config.filterNames = ('g', 'r', 'i', 'i2', 'z', 'y')
# Each filter maps onto a "standard filter".
# Here, both i and i2 map onto i2 which will be the standard
config.stdFilterNames = ('g', 'r', 'i2', 'i2', 'z', 'y')
# Pre-generated atmosphere table distributed with FGCM
config.atmosphereTableName = 'fgcm_atm_subaru1'

