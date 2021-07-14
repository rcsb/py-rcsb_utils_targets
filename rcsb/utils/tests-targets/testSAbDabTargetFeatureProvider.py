##
# File:    SAbDabTargetFeatureProviderTests.py
# Author:  J. Westbrook
# Date:    11-Jun-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing SAbDab target feature data.
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

from rcsb.utils.targets.SAbDabTargetFeatureProvider import SAbDabTargetFeatureProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class SAbDabTargetFeatureProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__seqMatchResultsPath = os.path.join(HERE, "test-data", "sabdab-vs-pdbprent-filtered-results.json.gz")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildSAbDabTargetsFeatures(self):
        stfP = SAbDabTargetFeatureProvider(cachePath=self.__cachePath, useCache=False)
        ok = stfP.testCache()
        self.assertFalse(ok)
        ok = stfP.buildFeatureList(self.__seqMatchResultsPath)
        self.assertTrue(ok)
        stfP = SAbDabTargetFeatureProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.testCache()
        self.assertTrue(ok)
        ok = stfP.hasFeatures("4g6k_2")
        self.assertTrue(ok)
        fL = stfP.getFeatures("4g6k_2")
        self.assertGreaterEqual(len(fL), 3)
        # "5bk0|A"
        ok = stfP.hasAssignment("5bk0.A")
        self.assertTrue(ok)
        ok = stfP.getAssignment("5bk0.A", "antigen_name")
        self.assertTrue(ok)


def buildSAbDabFeaturesTargets():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SAbDabTargetFeatureProviderTests("testBuildSAbDabTargetsFeatures"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = buildSAbDabFeaturesTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
