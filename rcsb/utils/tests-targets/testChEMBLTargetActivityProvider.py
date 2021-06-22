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

from rcsb.utils.targets.ChEMBLTargetActivityProvider import ChEMBLTargetActivityProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ChEMBLTargetActivityProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__mU = MarshalUtil(workPath=self.__cachePath)

    def tearDown(self):
        pass

    def testFetchMoleculeData(self):
        try:
            ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=False)
            ok = ctP.testCache()
            self.assertTrue(ok)
            chemblId = "CHEMBL1421"
            name, inchiKey, _ = ctP.getMoleculeDetails(chemblId)
            self.assertEqual(name, "DASATINIB")
            self.assertEqual(inchiKey, "ZBNZXTGUTAYRHI-UHFFFAOYSA-N")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchActivityDataMulti(self):
        try:
            ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=False)
            ok = ctP.testCache()
            self.assertTrue(ok)
            #
            tL = ["CHEMBL1987", "CHEMBL3243"]
            ok = ctP.fetchTargetActivityDataMulti(tL, numProc=2)
            self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchActivityData(self):
        try:
            ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=False)
            ok = ctP.testCache()
            self.assertTrue(ok)

            mL = ["CHEMBL200117"]
            for chemblId in mL:
                molD = ctP.getMoleculeDetails(chemblId)
                logger.info("%s molD %r", chemblId, molD)
            #
            tL = ["CHEMBL1987", "CHEMBL3243"]
            ok = ctP.fetchTargetActivityData(tL)
            self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def targetActivitySuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChEMBLTargetActivityProviderTests("testFetchMoleculeData"))
    suiteSelect.addTest(ChEMBLTargetActivityProviderTests("testFetchActivityData"))
    suiteSelect.addTest(ChEMBLTargetActivityProviderTests("testFetchActivityDataMulti"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = targetActivitySuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
