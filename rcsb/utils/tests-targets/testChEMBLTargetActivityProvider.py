##
# File:    ChEMBLTargetActivityProviderTests.py
# Author:  J. Westbrook
# Date:    15-Jun-2021
#
# Update:
#
#
##
"""
Tests for accessors for ChEMBL target activity data.
"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from chembl_webresource_client.settings import Settings

from rcsb.utils.targets.ChEMBLTargetActivityProvider import ChEMBLTargetActivityProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

Settings.Instance().TIMEOUT = 10  # pylint: disable=no-member
Settings.Instance().MAX_LIMIT = 50  # pylint: disable=no-member
Settings.MAX_LIMIT = 50


logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChEMBLTargetActivityProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__mU = MarshalUtil(workPath=self.__cachePath)

    def tearDown(self):
        pass

    def testFetchActivityData(self):
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=True)
            ok = ctP.testCache()
            self.assertTrue(ok)
            #
            tL = ["CHEMBL1987", "CHEMBL3243"]
            targetD = ctP.fetchActivityData(tL)
            logger.info("keys: %r", list(targetD.keys()))

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchMechanismData(self):
        oD = {}
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=True)
            ok = ctP.testCache()
            self.assertTrue(ok)
            tL = ["CHEMBL1987", "CHEMBL3243"]
            oD.update(ctP.getMechanismData(tL))
            #
            ok = self.__mU.doExport(os.path.join(self.__cachePath, "ChEMBL-targets", "chembl-target-mechanism.json"), oD, fmt="json", indent=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def targetActivitySuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChEMBLTargetActivityProviderTests("testFetchActivityData"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = targetActivitySuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
