##
# File:    SAbDabTargetProviderTests.py
# Author:  J. Westbrook
# Date:    30-Nov-2020
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing SAbDab target data.
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

from rcsb.utils.targets.SAbDabTargetProvider import SAbDabTargetProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class SAbDabTargetProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__fastaPath = os.path.join(HERE, "test-output", "sabdab-targets.fa")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testFetchSAbDabTargets(self):
        stP = SAbDabTargetProvider(cachePath=self.__cachePath, useCache=False)
        ok = stP.testCache()
        self.assertTrue(ok)
        ok = stP.exportFasta(self.__fastaPath)
        self.assertTrue(ok)


def fetchSAbDabTargets():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SAbDabTargetProviderTests("testFetchSAbDabTargets"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchSAbDabTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
