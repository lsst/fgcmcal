# See COPYRIGHT file at the top of the source tree.
#
# This file is part of fgcmcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Test the fgcmcal fgcmLoadReferenceCatalog code with testdata_jointcal/hsc.

"""

import unittest
import os
import numpy as np
import hpgeom as hpg
import esutil
import tempfile

import lsst.utils
import lsst.pipe.tasks
from lsst.meas.algorithms import ReferenceObjectLoader, LoadReferenceObjectsConfig

import lsst.fgcmcal as fgcmcal

import fgcmcalTestBase

ROOT = os.path.abspath(os.path.dirname(__file__))


class FgcmLoadReferenceTestHSC(fgcmcalTestBase.FgcmcalTestBase, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")
        try:
            lsst.utils.getPackageDir('obs_subaru')
        except LookupError:
            raise unittest.SkipTest("obs_subaru not setup")

        lsst.daf.butler.cli.cliLog.CliLog.initLog(longlog=False)

        cls.testDir = tempfile.mkdtemp(dir=ROOT, prefix="TestFgcm-")

        cls._importRepository('lsst.obs.subaru.HyperSuprimeCam',
                              os.path.join(cls.dataDir, 'hsc/repo'),
                              os.path.join(cls.dataDir, 'hsc', 'exports.yaml'))

    def test_fgcmLoadReference(self):
        """
        Test loading of the fgcm reference catalogs.
        """

        filterList = ['HSC-R', 'HSC-I']

        config = fgcmcal.FgcmLoadReferenceCatalogConfig()
        config.applyColorTerms = True
        config.filterMap = {'HSC-R': 'r', 'HSC-I': 'i'}
        config.colorterms.data = {}
        config.colorterms.data['ps1*'] = lsst.pipe.tasks.colorterms.ColortermDict()
        config.colorterms.data['ps1*'].data = {}
        config.colorterms.data['ps1*'].data['HSC-R'] = lsst.pipe.tasks.colorterms.Colorterm()
        config.colorterms.data['ps1*'].data['HSC-R'].primary = 'r'
        config.colorterms.data['ps1*'].data['HSC-R'].secondary = 'i'
        config.colorterms.data['ps1*'].data['HSC-R'].c0 = -0.000144
        config.colorterms.data['ps1*'].data['HSC-R'].c1 = 0.001369
        config.colorterms.data['ps1*'].data['HSC-R'].c2 = -0.008380
        config.colorterms.data['ps1*'].data['HSC-I'] = lsst.pipe.tasks.colorterms.Colorterm()
        config.colorterms.data['ps1*'].data['HSC-I'].primary = 'i'
        config.colorterms.data['ps1*'].data['HSC-I'].secondary = 'z'
        config.colorterms.data['ps1*'].data['HSC-I'].c0 = 0.000643
        config.colorterms.data['ps1*'].data['HSC-I'].c1 = -0.130078
        config.colorterms.data['ps1*'].data['HSC-I'].c2 = -0.006855

        refCatName = 'ps1_pv3_3pi_20170110'

        butler = lsst.daf.butler.Butler(self.repo, instrument='HSC', collections=['HSC/testdata',
                                                                                  'refcats/gen2'])
        refs = set(butler.registry.queryDatasets(refCatName))
        dataIds = [butler.registry.expandDataId(ref.dataId) for ref in refs]
        refCats = [butler.getDeferred(ref) for ref in refs]

        refConfig = LoadReferenceObjectsConfig()
        refConfig.filterMap = config.filterMap

        refObjLoader = ReferenceObjectLoader(dataIds=dataIds,
                                             refCats=refCats,
                                             name=refCatName,
                                             config=refConfig)

        loadCat = fgcmcal.FgcmLoadReferenceCatalogTask(refObjLoader=refObjLoader,
                                                       refCatName=refCatName,
                                                       config=config)

        ra = 337.656174
        dec = 0.823595
        rad = 0.1

        refCat = loadCat.getFgcmReferenceStarsSkyCircle(ra, dec, rad, filterList)

        # Check the number of mags and ranges
        self.assertEqual(len(filterList), refCat['refMag'].shape[1])
        self.assertEqual(len(filterList), refCat['refMagErr'].shape[1])
        self.assertLess(np.max(refCat['refMag'][:, 0]), 99.1)
        self.assertLess(np.max(refCat['refMag'][:, 1]), 99.1)
        self.assertLess(np.max(refCat['refMagErr'][:, 0]), 99.1)
        self.assertLess(np.max(refCat['refMagErr'][:, 1]), 99.1)
        test, = np.where((refCat['refMag'][:, 0] < 30.0)
                         & (refCat['refMag'][:, 1] < 30.0))
        self.assertGreater(test.size, 0)

        # Check the separations from the center
        self.assertLess(np.max(esutil.coords.sphdist(ra, dec, refCat['ra'], refCat['dec'])), rad)

        # And load a healpixel
        nside = 256
        pixel = 387520

        refCat = loadCat.getFgcmReferenceStarsHealpix(nside, pixel, filterList)

        ipring = hpg.angle_to_pixel(nside, refCat['ra'], refCat['dec'], nest=False)
        self.assertEqual(pixel, np.max(ipring))
        self.assertEqual(pixel, np.min(ipring))

    def test_fgcmLoadReferenceOtherFilters(self):
        """
        Test loading of the fgcm reference catalogs using unmatched filter names.
        """

        filterList = ['HSC-R2', 'HSC-I2']

        config = fgcmcal.FgcmLoadReferenceCatalogConfig()
        config.applyColorTerms = True
        config.filterMap = {'HSC-R2': 'r', 'HSC-I2': 'i'}
        config.colorterms.data = {}
        config.colorterms.data['ps1*'] = lsst.pipe.tasks.colorterms.ColortermDict()
        config.colorterms.data['ps1*'].data = {}
        config.colorterms.data['ps1*'].data['HSC-R2'] = lsst.pipe.tasks.colorterms.Colorterm()
        config.colorterms.data['ps1*'].data['HSC-R2'].primary = 'r'
        config.colorterms.data['ps1*'].data['HSC-R2'].secondary = 'i'
        config.colorterms.data['ps1*'].data['HSC-R2'].c0 = -0.000032
        config.colorterms.data['ps1*'].data['HSC-R2'].c1 = -0.002866
        config.colorterms.data['ps1*'].data['HSC-R2'].c2 = -0.012638
        config.colorterms.data['ps1*'].data['HSC-I2'] = lsst.pipe.tasks.colorterms.Colorterm()
        config.colorterms.data['ps1*'].data['HSC-I2'].primary = 'i'
        config.colorterms.data['ps1*'].data['HSC-I2'].secondary = 'z'
        config.colorterms.data['ps1*'].data['HSC-I2'].c0 = 0.001625
        config.colorterms.data['ps1*'].data['HSC-I2'].c1 = -0.200406
        config.colorterms.data['ps1*'].data['HSC-I2'].c2 = -0.013666

        refCatName = 'ps1_pv3_3pi_20170110'

        butler = lsst.daf.butler.Butler(self.repo, instrument='HSC', collections=['HSC/testdata',
                                                                                  'refcats/gen2'])
        refs = set(butler.registry.queryDatasets(refCatName))
        dataIds = [butler.registry.expandDataId(ref.dataId) for ref in refs]
        refCats = [butler.getDeferred(ref) for ref in refs]

        refConfig = LoadReferenceObjectsConfig()
        refConfig.filterMap = config.filterMap

        refObjLoader = ReferenceObjectLoader(dataIds=dataIds,
                                             refCats=refCats,
                                             name=refCatName,
                                             config=refConfig)

        loadCat = fgcmcal.FgcmLoadReferenceCatalogTask(refObjLoader=refObjLoader,
                                                       refCatName=refCatName,
                                                       config=config)

        ra = 337.656174
        dec = 0.823595
        rad = 0.1

        refCat = loadCat.getFgcmReferenceStarsSkyCircle(ra, dec, rad, filterList)

        self.assertEqual(len(filterList), refCat['refMag'].shape[1])
        self.assertEqual(len(filterList), refCat['refMagErr'].shape[1])
        test, = np.where((refCat['refMag'][:, 0] < 30.0)
                         & (refCat['refMag'][:, 1] < 30.0))
        self.assertGreater(test.size, 0)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
