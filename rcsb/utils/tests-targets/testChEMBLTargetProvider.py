##
# File:    ChEMBLTargetProviderTests.py
# Author:  J. Westbrook
# Date:    9-Nov-2020
#
# Update:
#
#
##
"""
Tests for accessors for ChEMBL target data.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.targets.ChEMBLTargetProvider import ChEMBLTargetProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChEMBLTargetProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mU = MarshalUtil(workPath=self.__cachePath)

    def tearDown(self):
        pass

    def testFetchChEMBLTargets(self):
        try:
            pcP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=False)
            ok = pcP.testCache()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def chemblTargetSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChEMBLTargetProviderTests("testFetchChEMBLTargets"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = chemblTargetSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
