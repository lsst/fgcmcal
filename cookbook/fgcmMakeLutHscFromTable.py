"""
HSC-specific overrides for FgcmMakeLut
"""

# Short-code Filter names
config.filterNames = ('g', 'r', 'i', 'z', 'y')
# Each filter maps onto a "standard filter".
config.stdFilterNames = ('g', 'r', 'i', 'z', 'y')
# Pre-generated atmosphere table distributed with FGCM
config.atmosphereTableName = 'fgcm_atm_subaru3'

