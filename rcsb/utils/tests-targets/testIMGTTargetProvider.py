##
# File:    IMGTTargetProviderTests.py
# Author:  J. Westbrook
# Date:    3-Jul-2021
# Version: 0.001
#
# Update:
#  26-Oct-2026 dwp Add test for IMGT source entry IDs (fails once IMGT uses extended PDB IDs)
#
#
##
"""
Tests for processing IMGT reference sequence data.
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

from rcsb.utils.targets.IMGTTargetProvider import IMGTTargetProvider

# from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


def getNonShortEntryIds(entryIds):
    """Return the entry IDs that aren't short (4-character) PDB IDs, e.g. extended (12-character) PDB IDs."""
    return sorted({eId for eId in entryIds if len(eId) != 4})


class IMGTTargetProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    @unittest.skipIf(True, "Avoid re-downloading large file - should have already been downloaded by testIMGTTargetFeatureProvider.py")
    def testAPrepTargetData(self):
        imgtP = IMGTTargetProvider(self.__cachePath, useCache=False)
        ok = imgtP.testCache()
        self.assertTrue(ok)

    def testBExportFasta(self):
        imgtP = IMGTTargetProvider(self.__cachePath, useCache=True)
        ok = imgtP.testCache()
        self.assertTrue(ok)
        ok = imgtP.exportFasta()
        self.assertTrue(ok)

    def testCSourceEntryIds(self):
        """Check that the IMGT source data (REMARK 410 chain IDs, '<entryId>_<authAsymId>') still uses short (4-character) PDB IDs.

        Like testBExportFasta, reads the cache built from the source data by testIMGTTargetFeatureProvider.py.
        """
        imgtP = IMGTTargetProvider(self.__cachePath, useCache=True)
        self.assertTrue(imgtP.testCache())
        # Fails once the source data uses extended PDB IDs; then update this package and rcsb.utils.dictionary for them
        badL = getNonShortEntryIds({"_".join(chainId.split("_")[:-1]) for chainId in imgtP.getChains()})
        self.assertEqual(badL, [], "IMGT source data has entry IDs that aren't 4-char PDB IDs: %r" % badL[:10])


def imgtTargetPrep():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(IMGTTargetProviderTests("testAPrepTargetData"))
    suiteSelect.addTest(IMGTTargetProviderTests("testBExportFasta"))
    suiteSelect.addTest(IMGTTargetProviderTests("testCSourceEntryIds"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = imgtTargetPrep()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
