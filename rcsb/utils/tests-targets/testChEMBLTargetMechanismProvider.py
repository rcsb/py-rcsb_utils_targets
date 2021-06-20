##
# File:    ChEMBLTargetMechanismProviderTests.py
# Author:  J. Westbrook
# Date:    15-Jun-2021
#
# Update:
#
#
##
"""
Tests for accessors for ChEMBL target mechanism data.
"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.targets.ChEMBLTargetMechanismProvider import ChEMBLTargetMechanismProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChEMBLTargetMechanismProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__mU = MarshalUtil(workPath=self.__cachePath)

    def tearDown(self):
        pass

    def testFetchMechanismData(self):
        try:
            ctP = ChEMBLTargetMechanismProvider(cachePath=self.__cachePath, useCache=True)
            ok = ctP.testCache()
            self.assertTrue(ok)
            #
            tL = ["CHEMBL1987", "CHEMBL3243"]
            targetD = ctP.fetchMechanismData(tL)
            logger.info("keys: %r", list(targetD.keys()))

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def targetMechanismSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChEMBLTargetMechanismProviderTests("testFetchMechanismData"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = targetMechanismSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
