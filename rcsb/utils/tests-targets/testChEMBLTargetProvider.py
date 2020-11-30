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
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mU = MarshalUtil(workPath=self.__cachePath)

    def tearDown(self):
        pass

    def testGetUniChemSources(self):
        try:
            srL = self.__mU.doImport(os.path.join(self.__dataPath, "UC_SOURCE.tdd"), fmt="tdd", rowFormat="list")
            # logger.info("srDL %r", srL)
            uD = {}
            for sr in srL:
                if len(sr) < 12:
                    logger.info("bad sr %r", sr)
                    continue
                uD[int(sr[0])] = {"name": sr[1], "baseUrl": sr[10], "entryUrl": sr[11]}
            logger.info("uD = %r", uD)
            ok = self.__mU.doExport(os.path.join(self.__cachePath, "unichem-source-list.json"), uD, fmt="json", indent=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchChEMBLTargets(self):
        try:
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=False, addTaxonomy=False)
            ok = ctP.testCache()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testActivityData(self):
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=True, addTaxonomy=False)
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

    def testMechanismData(self):
        oD = {}
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=True, addTaxonomy=False)
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

    def testDrugData(self):
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=True, addTaxonomy=False)
            ok = ctP.testCache()
            self.assertTrue(ok)
            mL = ["CHEMBL815", "CHEMBL426559", "CHEMBL548", "CHEMBL2442750", "CHEMBL1201379", "CHEMBL3039498", "CHEMBL3137332"]
            #
            oD = ctP.getDrugData(mL)
            #
            logger.info("molD %r", list(oD.keys()))
            ok = self.__mU.doExport(os.path.join(self.__cachePath, "ChEMBL-targets", "chembl-drug.json"), oD, fmt="json", indent=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMoleculeData(self):
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=True, addTaxonomy=False)
            ok = ctP.testCache()
            self.assertTrue(ok)
            mL = ["CHEMBL815", "CHEMBL426559", "CHEMBL548", "CHEMBL2442750", "CHEMBL1201379", "CHEMBL3039498", "CHEMBL3137332"]
            #
            oD = ctP.getMoleculeData(mL)
            #
            logger.info("molD %r", list(oD.keys()))
            ok = self.__mU.doExport(os.path.join(self.__cachePath, "ChEMBL-targets", "chembl-molecule.json"), oD, fmt="json", indent=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testMolelcuesByInChI(self):
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=True, addTaxonomy=False)
            ok = ctP.testCache()
            self.assertTrue(ok)
            mL = ["NXFFJDQHYLNEJK-CYBMUJFWSA-N", "WWSWYXNVCBLWNZ-QIZQQNKQSA-N"]
            oD = ctP.getMoleculeByInChIKey(mL)
            #
            ok = self.__mU.doExport(os.path.join(self.__cachePath, "ChEMBL-targets", "chembl-inchikey-matches.json"), oD, fmt="json", indent=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testUniChemData(self):
        try:
            logger.info("MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=True, addTaxonomy=False)
            ok = ctP.testCache()
            self.assertTrue(ok)
            mL = ["NXFFJDQHYLNEJK-CYBMUJFWSA-N", "WWSWYXNVCBLWNZ-QIZQQNKQSA-N"]
            oD = ctP.getUniChemData(mL)
            #
            ok = self.__mU.doExport(os.path.join(self.__cachePath, "ChEMBL-targets", "unichem-inchikey-matches.json"), oD, fmt="json", indent=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    #
    # @unittest.skipIf(skipFull, "Very long test")
    def testFetchChEMBLTargetsWithTax(self):
        try:
            ctP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=False, addTaxonomy=True)
            ok = ctP.testCache()
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
