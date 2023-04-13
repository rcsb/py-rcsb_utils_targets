##
# File:    CARDTargetAnnotationProviderTests.py
# Author:  D. Piehl
# Date:    6-Mar-2023
# Version: 0.001
#
# Update:
#
##
"""
Tests for utilities managing CARD target annotation data.
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

from rcsb.utils.targets.CARDTargetAnnotationProvider import CARDTargetAnnotationProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class CARDTargetAnnotationProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__seqMatchResultsPath = os.path.join(HERE, "test-data", "card-vs-pdbprent-filtered-results.json.gz")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildCARDTargetsAnnotations(self):
        stfP = CARDTargetAnnotationProvider(cachePath=self.__cachePath, useCache=False)
        ok = stfP.testCache()
        self.assertFalse(ok)
        ok = stfP.buildAnnotationList(self.__seqMatchResultsPath, useTaxonomy=False)
        # ok = stfP.buildAnnotationList(self.__seqMatchResultsPath, useTaxonomy=True)
        self.assertTrue(ok)
        stfP = CARDTargetAnnotationProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.testCache()
        self.assertTrue(ok)
        ok = stfP.hasAnnotation("5f64_1")
        self.assertTrue(ok)
        aD = stfP.getAnnotation("5f64_1")
        self.assertGreaterEqual(len(list(aD)), 10)
        tnL = stfP.getTreeNodeList()
        self.assertGreaterEqual(len(tnL), 100)


def buildCARDAnnotationsTargets():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CARDTargetAnnotationProviderTests("testBuildCARDTargetsAnnotations"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = buildCARDAnnotationsTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
