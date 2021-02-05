"""
HSC-specific overrides for FgcmMakeLut
"""
# Physical filter labels to generate LUT
config.physicalFilters = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y']
# For RC2, we must explicitly map HSC-R to HSC-R and HSC-I to HSC-I
# in the override configuration, as the obs_subaru default assumes both
# HSC-R + HSC-R2 and HSC-I + HSC-I2 data
config.stdPhysicalFilterOverrideMap = {'HSC-R': 'HSC-R',
                                       'HSC-I': 'HSC-I'}
# Pre-generated atmosphere table distributed with FGCM
config.atmosphereTableName = 'fgcm_atm_subaru3'

