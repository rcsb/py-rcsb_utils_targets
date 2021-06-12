##
# File:    CARDTargetFeatureProviderTests.py
# Author:  J. Westbrook
# Date:    11-Jun-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing CARD target data.
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

from rcsb.utils.targets.CARDTargetFeatureProvider import CARDTargetFeatureProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class CARDTargetProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__seqMatchResultsPath = os.path.join(HERE, "test-data", "pdb-card-results.json")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildCARDTargetsFeatures(self):
        stfP = CARDTargetFeatureProvider(cachePath=self.__cachePath, useCache=False)
        ok = stfP.testCache()
        self.assertFalse(ok)
        ok = stfP.buildFeatureList(self.__seqMatchResultsPath)
        self.assertTrue(ok)
        stfP = CARDTargetFeatureProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.testCache()
        self.assertTrue(ok)
        ok = stfP.hasFeatures("3g7e_1")
        self.assertTrue(ok)
        fL = stfP.getFeatures("3g7e_1")
        self.assertGreaterEqual(len(fL), 2)


def buildCARDFeaturesTargets():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CARDTargetProviderTests("testBuildCARDTargetsFeatures"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = buildCARDFeaturesTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
