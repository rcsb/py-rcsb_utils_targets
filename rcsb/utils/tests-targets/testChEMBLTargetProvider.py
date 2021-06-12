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

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from chembl_webresource_client.settings import Settings

from rcsb.utils.targets.ChEMBLTargetProvider import ChEMBLTargetProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

Settings.Instance().TIMEOUT = 10  # pylint: disable=no-member
Settings.Instance().MAX_LIMIT = 50  # pylint: disable=no-member
Settings.MAX_LIMIT = 50


logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChEMBLTargetProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__fastaPath = os.path.join(HERE, "test-output", "chembl-targets.fa")
        self.__taxonPath = os.path.join(HERE, "test-output", "chembl-targets-taxon.tdd")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mU = MarshalUtil(workPath=self.__cachePath)

    def tearDown(self):
        pass

    def testFetchChEMBLTargets(self):
        try:
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=False)
            ok = ctP.testCache()
            self.assertTrue(ok)
            ok = ctP.exportFasta(self.__fastaPath, self.__taxonPath, addTaxonomy=False)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchActivityData(self):
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=True)
            ok = ctP.testCache()
            self.assertTrue(ok)
            # P43088|CHEMBL1987|9606
            # P43093|CHEMBL6179|294748
            tL = ["CHEMBL1987", "CHEMBL6179"]
            targetD, molD = ctP.getActivityData(tL)
            logger.info("keys: %r", list(molD.keys()))
            ok = self.__mU.doExport(os.path.join(self.__cachePath, "ChEMBL-targets", "chembl-target-activity.json"), targetD, fmt="json", indent=3)
            ok = self.__mU.doExport(os.path.join(self.__cachePath, "ChEMBL-targets", "chembl-molecule-activity.json"), molD, fmt="json", indent=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchMechanismData(self):
        oD = {}
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=True)
            ok = ctP.testCache()
            self.assertTrue(ok)
            # P43088|CHEMBL1987|9606
            # P43093|CHEMBL6179|294748
            tL = ["CHEMBL1987", "CHEMBL6179"]
            oD.update(ctP.getMechanismData(tL))
            #
            ok = self.__mU.doExport(os.path.join(self.__cachePath, "ChEMBL-targets", "chembl-target-mechanism.json"), oD, fmt="json", indent=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    #
    @unittest.skipIf(skipFull, "Very long test")
    def testFetchChEMBLTargetsWithTax(self):
        try:
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=True)
            ok = ctP.testCache()
            self.assertTrue(ok)
            ok = ctP.exportFasta(self.__fastaPath, self.__taxonPath, addTaxonomy=True)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def chemblTargetSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChEMBLTargetProviderTests("testFetchChEMBLTargets"))
    suiteSelect.addTest(ChEMBLTargetProviderTests("testFetchMechanismData"))
    suiteSelect.addTest(ChEMBLTargetProviderTests("testFetchActivityData"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = chemblTargetSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
