config.connections.cycleNumber = 2

# Flat metadata is not available in testdata_jointcal
config.use_flat_metadata = False
config.approximate_wcs_jacobian = False
config.epoch_time =  "2023-08-22"
config.physical_filters = [
    "SDSSg_65mm~empty",
    "SDSSr_65mm~empty",
    "SDSSi_65mm~empty",
]
