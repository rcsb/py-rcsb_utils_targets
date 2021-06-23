##
# File:    DrugBankTargetProviderTests.py
# Author:  J. Westbrook
# Date:    7-Nov-2020
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing DrugBank target data.

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

from rcsb.utils.targets.DrugBankTargetCofactorProvider import DrugBankTargetCofactorProvider
from rcsb.utils.targets.DrugBankTargetProvider import DrugBankTargetProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class DrugBankTargetProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        configPath = os.path.join(HERE, "test-data", "drugbank-config-example.yml")
        configName = "site_info_configuration"
        cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
        self.__user = cfgOb.get("_DRUGBANK_AUTH_USERNAME", sectionName=configName)
        self.__pw = cfgOb.get("_DRUGBANK_AUTH_PASSWORD", sectionName=configName)
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__fastaPath = os.path.join(HERE, "test-output", "drugbank-targets.fa")
        self.__taxonPath = os.path.join(HERE, "test-output", "drugbank-targets-taxon.tdd")
        #
        self.__seqMatchResultsPath = os.path.join(HERE, "test-data", "drugbank-vs-pdbprent-filtered-results.json.gz")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    @unittest.skipIf(skipFull, "Very long test")
    def testAAFetchDrugBankTargets(self):
        dbtP = DrugBankTargetProvider(cachePath=self.__cachePath, useCache=False, username=self.__user, password=self.__pw)
        ok = dbtP.testCache()
        self.assertTrue(ok)
        ok = dbtP.exportFasta(self.__fastaPath, self.__taxonPath, addTaxonomy=False)
        self.assertTrue(ok)

    @unittest.skipIf(skipFull, "Very long test")
    def testBBFetchDrugBankTargetsWithTaxonomy(self):
        dbtP = DrugBankTargetProvider(cachePath=self.__cachePath, useCache=False, username=self.__user, password=self.__pw, addTaxonomy=True)
        ok = dbtP.testCache()
        self.assertTrue(ok)
        ok = dbtP.exportFasta(self.__fastaPath, self.__taxonPath, addTaxonomy=True)
        self.assertTrue(ok)

    @unittest.skipIf(skipFull, "Very long test")
    def testQQBuildDrugBankTargetsFeatures(self):
        #
        stfP = DrugBankTargetCofactorProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.buildCofactorList(self.__seqMatchResultsPath)
        self.assertTrue(ok)
        stfP = DrugBankTargetCofactorProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.testCache()
        self.assertTrue(ok)
        ok = stfP.hasCofactor("5vbn_2")
        self.assertTrue(ok)
        fL = stfP.getCofactors("5vbn_2")
        self.assertGreaterEqual(len(fL), 1)


def fetchDrugBankTargets():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DrugBankTargetProviderTests("testAAFetchDrugBankTargets"))
    suiteSelect.addTest(DrugBankTargetProviderTests("testBBFetchDrugBankTargetsWithTaxonomy"))
    suiteSelect.addTest(DrugBankTargetProviderTests("testQQBuildDrugBankTargetsFeatures"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchDrugBankTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
