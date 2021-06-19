##
# File:    ChEMBLTargetCofactorProviderTests.py
# Author:  J. Westbrook
# Date:    15-Jun-2021
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
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.targets.ChEMBLTargetCofactorProvider import ChEMBLTargetCofactorProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ChEMBLTargetCofactorProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__seqMatchResultsPath = os.path.join(HERE, "test-data", "chembl-vs-pdbprent-filtered-results.json.gz")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildChEMBLTargetsCofactors(self):
        stfP = ChEMBLTargetCofactorProvider(cachePath=self.__cachePath, useCache=False)
        ok = stfP.testCache()
        self.assertFalse(ok)
        ok = stfP.buildCofactorList(self.__seqMatchResultsPath)
        self.assertTrue(ok)
        stfP = ChEMBLTargetCofactorProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.testCache()
        self.assertTrue(ok)
        ok = stfP.hasCofactor("5fn7_1")
        self.assertTrue(ok)
        aL = stfP.getCofactors("5fn7_1")
        self.assertGreaterEqual(len(aL[0]["cofactors"]), 5)


def buildChEMBLTargetCofactors():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChEMBLTargetCofactorProviderTests("testBuildChEMBLTargetsCofactors"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = buildChEMBLTargetCofactors()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
