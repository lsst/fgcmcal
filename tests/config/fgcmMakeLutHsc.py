config.physicalFilters = ['HSC-G', 'HSC-R', 'HSC-I']
# We need the physical filter overrides here to override obs_subaru configs
config.stdPhysicalFilterOverrideMap = {'HSC-R': 'HSC-R',
                                       'HSC-I': 'HSC-I'}
config.atmosphereTableName = 'fgcm_atm_subaru2_test'
