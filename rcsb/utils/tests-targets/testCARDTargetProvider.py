##
# File:    CARDTargetProviderTests.py
# Author:  J. Westbrook
# Date:    27-Nov-2020
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

from rcsb.utils.targets.CARDTargetProvider import CARDTargetProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class CARDTargetProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__fastaPath = os.path.join(HERE, "test-output", "card-targets.fa")
        self.__taxonPath = os.path.join(HERE, "test-output", "card-targets-taxon.tdd")
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testFetchCARDTargets(self):
        ctP = CARDTargetProvider(cachePath=self.__cachePath, useCache=False)
        ok = ctP.testCache()
        self.assertTrue(ok)
        ok = ctP.exportCardFasta(self.__fastaPath, self.__taxonPath)
        self.assertTrue(ok)
        aroIdL = ["ARO:3000096", "ARO:3000015"]
        for aroId in aroIdL:
            lineageL = ctP.getLineage(aroId)
            logger.info("Lineage list for aroId %s: %r", aroId, lineageL)
            self.assertGreaterEqual(len(lineageL), 2)


def fetchCARDTargets():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CARDTargetProviderTests("testFetchCARDTargets"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchCARDTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
