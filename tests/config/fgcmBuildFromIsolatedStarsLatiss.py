import os
from lsst.obs.lsst.filters import LATISS_FILTER_DEFINITIONS

# The filterMap and bands are for the small subset of bands used in the tests
config.physicalFilterMap = {
    "SDSSg_65mm~empty": "g",
    "SDSSr_65mm~empty": "r",
    "SDSSi_65mm~empty": "i",
}
config.requiredBands = ["g", "r"]
config.primaryBands = ["r"]
# The coarseNside is set appropriate to the area of the test data
config.coarseNside = 64
# The tests are done with only the brightest reference stars
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.minimum = 50.0

config.instFluxField = "apFlux_35_0_instFlux"
config.apertureInnerInstFluxField = "apFlux_35_0_instFlux"
config.apertureOuterInstFluxField = "apFlux_50_0_instFlux"

config.minPerBand = 2
config.connections.ref_cat = "atlas_refcat2_20220201"

configDir = os.path.join(os.path.dirname(__file__))
config.physicalFilterMap = LATISS_FILTER_DEFINITIONS.physical_to_band
config.doSubtractLocalBackground = True
config.sourceSelector["science"].flags.bad.append("localBackground_flag")
config.sourceSelector["science"].signalToNoise.fluxField = "apFlux_35_0_instFlux"
config.sourceSelector["science"].signalToNoise.errField = "apFlux_35_0_instFluxErr"
config.sourceSelector["science"].signalToNoise.minimum = 11.0
config.fgcmLoadReferenceCatalog.load(os.path.join(configDir, "filterMapLatiss.py"))
config.fgcmLoadReferenceCatalog.applyColorTerms = True
config.fgcmLoadReferenceCatalog.colorterms.load(os.path.join(configDir, "colortermsLatiss.py"))
config.fgcmLoadReferenceCatalog.referenceSelector.doSignalToNoise = True
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.fluxField = "i_flux"
config.fgcmLoadReferenceCatalog.referenceSelector.signalToNoise.errField = "i_fluxErr"
