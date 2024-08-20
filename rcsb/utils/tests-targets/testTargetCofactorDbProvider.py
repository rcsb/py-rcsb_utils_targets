##
# File:    testTargetCofactorDbProvider.py
# Author:  Dennis Piehl
# Date:    20-Aug-2024
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing ChEMBL target cofactor data.
"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.targets.PharosTargetCofactorProvider import PharosTargetCofactorProvider
from rcsb.utils.targets.TargetCofactorDbProvider import TargetCofactorDbProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class TargetCofactorDbProviderTests(unittest.TestCase):
    skipTests = False

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__configPath = os.path.join(self.__mockTopPath, "config", "dbload-setup-example.yml")
        self.__configName = "site_info_configuration"
        self.__cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=self.__configName, mockTopPath=self.__mockTopPath)
        #
        self.__seqMatchResultsPath = os.path.join(HERE, "test-data", "pharos-vs-pdbprent-filtered-results.json.gz")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    @unittest.skipIf(skipTests, "Partly redundant with testPharosTargetCofactorProvider")
    def testTargetsCofactorsDbPharos(self):
        resourceName = "pharos"
        #
        # First re-load (or re-build) the resource
        stfP = PharosTargetCofactorProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.testCache()
        if not ok:
            logger.info("Resource cache not loadable--rebuilding")
            ok = stfP.buildCofactorList(self.__seqMatchResultsPath)
            self.assertTrue(ok)
            stfP = PharosTargetCofactorProvider(cachePath=self.__cachePath, useCache=True)
            ok = stfP.testCache()
        self.assertTrue(ok)
        ok = stfP.hasTarget("5fn7_1")
        self.assertTrue(ok)
        aL = stfP.getTargets("5fn7_1")
        self.assertGreaterEqual(len(aL[0]["cofactors"]), 5)
        #
        # Next test the loading of data to mongo
        tcDbP = TargetCofactorDbProvider(cachePath=self.__cachePath, cfgOb=self.__cfgOb, cofactorResourceName=resourceName)
        okLoad = tcDbP.loadCofactorData(resourceName, stfP)
        logger.info("%r cofactor data DB load status (%r)", resourceName, okLoad)
        numLoaded = tcDbP.cofactorDbCount()
        logger.info("%r cofactor data DB load count (%r)", resourceName, numLoaded)
        ok = numLoaded > 5
        self.assertTrue(ok)
        #
        # Now test access of data in mongo
        tDL = tcDbP.fetchCofactorData("5fn7_1")
        ok = len(tDL) > 0
        self.assertTrue(ok)


def testTargetCofactorDbProviderSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TargetCofactorDbProviderTests("testTargetsCofactorsDbPharos"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = testTargetCofactorDbProviderSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
