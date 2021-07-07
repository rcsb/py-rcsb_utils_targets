##
# File:    IMGTTargetFeatureProviderTests.py
# Author:  J. Westbrook
# Date:    11-Jun-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing IMGT target feature data.
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

from rcsb.utils.targets.IMGTTargetFeatureProvider import IMGTTargetFeatureProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class IMGTTargetFeatureProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildIMGTTargetsFeatures(self):
        imgtP = IMGTTargetFeatureProvider(cachePath=self.__cachePath, useCache=False)
        ok = imgtP.testCache()
        self.assertFalse(ok)
        ok = imgtP.buildFeatureList()
        self.assertTrue(ok)
        imgtP = IMGTTargetFeatureProvider(cachePath=self.__cachePath, useCache=True)
        ok = imgtP.testCache()
        self.assertTrue(ok)
        ok = imgtP.hasFeatures("6bla.L")
        self.assertTrue(ok)
        fL = imgtP.getFeatures("6bla.L")
        self.assertGreaterEqual(len(fL), 11)
        logger.debug("(%d) %r", len(fL), fL)


def buildIMGTFeaturesTargets():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(IMGTTargetFeatureProviderTests("testBuildIMGTTargetsFeatures"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = buildIMGTFeaturesTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
