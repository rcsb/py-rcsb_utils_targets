##
# File:    PharosTargetProviderTests.py
# Author:  J. Westbrook
# Date:    9-Nov-2020
#
# Update:
#
#
##
"""
Tests for accessors for Pharos target data.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.targets.PharosTargetProvider import PharosTargetProvider
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PharosTargetProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        configPath = os.path.join(HERE, "test-data", "pharos-config-example.yml")
        configName = "site_info_configuration"
        cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
        self.__user = cfgOb.get("_MYSQL_DB_USER_NAME", sectionName=configName)
        self.__pw = cfgOb.get("_MYSQL_DB_PASSWORD", sectionName=configName)
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dirPath = os.path.join(self.__cachePath, "Pharos-targets")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        self.__pharosFixture()

    def tearDown(self):
        pass

    def __pharosFixture(self):
        try:
            ok = False
            fU = FileUtil()
            srcPath = os.path.join(self.__dataPath, "Pharos")
            dstPath = self.__dirPath
            for fn in ["drug_activity", "cmpd_activity", "target", "protein", "t2tc"]:
                inpPath = os.path.join(srcPath, fn + ".tdd.gz")
                outPath = os.path.join(dstPath, fn + ".tdd.gz")
                fU.get(inpPath, outPath)
                fU.uncompress(outPath, outputDir=dstPath)
            ok = True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        return ok

    @unittest.skipIf(skipFull, "Very long test")
    def testFetchAndLoadPharosTargets(self):
        try:
            ptP = PharosTargetProvider(cachePath=self.__cachePath, useCache=False, reloadDb=True, mysqlUser=self.__user, mysqlPassword=self.__pw)
            ok = ptP.testCache()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFull, "Very long test")
    def testExportPharosTargets(self):
        try:
            ptP = PharosTargetProvider(cachePath=self.__cachePath, useCache=True, reloadDb=False, mysqlUser=self.__user, mysqlPassword=self.__pw)
            ok = ptP.testCache()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testExportPharosTargetFasta(self):
        try:
            ptP = PharosTargetProvider(cachePath=self.__cachePath, useCache=True, reloadDb=False)
            ok = ptP.testCache()
            self.assertTrue(ok)
            fastaPath = os.path.join(self.__dirPath, "pharos_targets_notax.fa")
            ok = ptP.exportProteinFasta(fastaPath, addTaxonomy=False)
            self.assertTrue(ok)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFull, "Very long test")
    def testExportPharosTargetFastaTax(self):
        try:
            ptP = PharosTargetProvider(cachePath=self.__cachePath, useCache=True, reloadDb=False)
            ok = ptP.testCache()
            self.assertTrue(ok)
            #
            fastaPath = os.path.join(self.__dirPath, "pharos_targets_tax.fa")
            ok = ptP.exportProteinFasta(fastaPath, addTaxonomy=True)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testExportPharosCofactors(self):
        try:
            ptP = PharosTargetProvider(cachePath=self.__cachePath, useCache=True, reloadDb=False)
            ok = ptP.testCache()
            self.assertTrue(ok)
            cofactorDataPath = os.path.join(self.__dirPath, "pharos_cofactors.json")
            ok = ptP.exportCofactors(cofactorDataPath, fmt="json")
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def pharosTargetSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PharosTargetProviderTests("testFetchPharosTargets"))
    suiteSelect.addTest(PharosTargetProviderTests("testFetchAndLoadPharosTargets"))
    suiteSelect.addTest(PharosTargetProviderTests("testExportPharosTargetFasta"))
    suiteSelect.addTest(PharosTargetProviderTests("testExportPharosCofactors"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = pharosTargetSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
