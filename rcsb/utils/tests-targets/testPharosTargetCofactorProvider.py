##
# File:    testPharosTargetCofactorProvider.py
# Author:  J. Westbrook
# Date:    15-Jun-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing Pharos target cofactor data.
"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.targets.PharosTargetCofactorProvider import PharosTargetCofactorProvider
from rcsb.utils.targets.PharosTargetCofactorProvider import PharosTargetCofactorAccessor

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PharosTargetCofactorProviderTests(unittest.TestCase):
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

    def testBuildPharosTargetsCofactors(self):
        stfP = PharosTargetCofactorProvider(cachePath=self.__cachePath, useCache=False)
        ok = stfP.testCache()
        self.assertFalse(ok)
        ok = stfP.buildCofactorList(self.__seqMatchResultsPath)
        self.assertTrue(ok)
        stfP = PharosTargetCofactorProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.testCache()
        self.assertTrue(ok)
        ok = stfP.hasTarget("5fn7_1")
        self.assertTrue(ok)
        aL = stfP.getTargets("5fn7_1")
        self.assertGreaterEqual(len(aL[0]["cofactors"]), 5)

    def testPharosTargetsCofactorsDb(self):
        # First test the loading of data to mongo
        stfP = PharosTargetCofactorProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.testCache()
        self.assertTrue(ok)
        ok = stfP.hasTarget("5fn7_1")
        self.assertTrue(ok)
        ok = stfP.loadCofactorData(cfgOb=self.__cfgOb)
        self.assertTrue(ok)
        #
        # Now test access of data in mongo
        stfA = PharosTargetCofactorAccessor(cachePath=self.__cachePath, cfgOb=self.__cfgOb)
        ok = stfA.testCache()
        tDL = stfA.getTargets("5fn7_1")
        ok = len(tDL) > 0
        self.assertTrue(ok)


def buildPharosTargetCofactors():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PharosTargetCofactorProviderTests("testBuildPharosTargetsCofactors"))
    suiteSelect.addTest(PharosTargetCofactorProviderTests("testPharosTargetsCofactorsDb"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = buildPharosTargetCofactors()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
