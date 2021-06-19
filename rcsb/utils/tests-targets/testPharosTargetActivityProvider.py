##
# File:    PharosTargetActivityProviderTests.py
# Author:  J. Westbrook
# Date:    17-Jun-2021
#
# Update:
#
#
##
"""
Tests for accessors for Pharos target activity data.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.targets.PharosTargetActivityProvider import PharosTargetActivityProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PharosTargetActivityProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dirPath = os.path.join(self.__cachePath, "Pharos-target-activity")
        self.__dataPath = os.path.join(HERE, "test-data")
        #
        self.__pharosFixture()

    def tearDown(self):
        pass
        #

    def __pharosFixture(self):
        try:
            ok = False
            fU = FileUtil()
            srcPath = os.path.join(self.__dataPath, "Pharos")
            dstPath = self.__dirPath
            for fn in ["drug_activity", "cmpd_activity"]:
                inpPath = os.path.join(srcPath, fn + ".tdd.gz")
                outPath = os.path.join(dstPath, fn + ".tdd.gz")
                fU.get(inpPath, outPath)
                fU.uncompress(outPath, outputDir=dstPath)
                fU.remove(outPath)
            ok = True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        return ok

    def testExportPharosActivity(self):
        try:
            ptP = PharosTargetActivityProvider(cachePath=self.__cachePath, useCache=True)
            ok = ptP.testCache()
            self.assertTrue(ok)
            ok = ptP.fetchTargetActivityData()
            self.assertTrue(ok)
            #
            ptP = PharosTargetActivityProvider(cachePath=self.__cachePath, useCache=True)
            ok = ptP.testCache(minCount=1)
            #
            ok = ptP.hasTargetActivity("3707")
            self.assertTrue(ok)
            aL = ptP.getTargetActivity("3707")
            self.assertGreater(len(aL), 125)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def pharosTargetActivitySuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PharosTargetActivityProviderTests("testExportPharosActivity"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = pharosTargetActivitySuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
