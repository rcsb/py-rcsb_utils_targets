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

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest
import requests

from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.targets.PharosTargetProvider import PharosTargetProvider
from rcsb.utils.io.FileUtil import FileUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PharosTargetProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        configPath = os.path.join(HERE, "test-data", "pharos-config-example.yml")
        self.__configName = "site_info_configuration"
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName)
        self.__user = self.__cfgOb.get("_MYSQL_DB_USER_NAME", sectionName=self.__configName)
        self.__pw = self.__cfgOb.get("_MYSQL_DB_PASSWORD", sectionName=self.__configName)
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dirPath = os.path.join(self.__cachePath, "Pharos-targets")
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
            for fn in ["drug_activity", "cmpd_activity", "target", "protein", "t2tc"]:
                inpPath = os.path.join(srcPath, fn + ".tdd.gz")
                outPath = os.path.join(dstPath, fn + ".tdd.gz")
                fU.get(inpPath, outPath)
                fU.uncompress(outPath, outputDir=dstPath)
                fU.remove(outPath)
            fU.put(os.path.join(srcPath, "pharos-readme.txt"), os.path.join(dstPath, "pharos-readme.txt"))
            ok = True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            ok = False
        return ok

    @unittest.skip("Bootstrap test")
    def testBootstrap(self):
        try:
            ptP = PharosTargetProvider(cachePath=self.__cachePath, useCache=False, reloadDb=False)
            configPath = os.path.join(TOPDIR, "rcsb", "mock-data", "config", "dbload-setup-example.yml")
            configName = "site_info_remote_configuration"
            cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
            ok = ptP.backup(cfgOb, configName)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFull, "Database dependency")
    def testFetchAndLoadPharosTargets(self):
        try:
            # Now about 630s on macos
            ptP = PharosTargetProvider(cachePath=self.__cachePath, useCache=False, reloadDb=True, fromDb=True, mysqlUser=self.__user, mysqlPassword=self.__pw)
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
            fastaPath = self.__cachePath = os.path.join(HERE, "test-output", "pharos-targets.fa")
            taxonPath = self.__cachePath = os.path.join(HERE, "test-output", "pharos-targets-taxon.tdd")
            ok = ptP.exportProteinFasta(fastaPath, taxonPath, addTaxonomy=False)
            self.assertTrue(ok)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFull, "Internal test")
    def testStashDependencies(self):
        try:
            ptP = PharosTargetProvider(cachePath=self.__cachePath, useCache=True, reloadDb=False, fromDb=False)
            ok = ptP.testCache()
            self.assertTrue(ok)
            #
            ok = ptP.backup(self.__cfgOb, self.__configName)
            self.assertTrue(ok)
            #
            ok = ptP.restore(self.__cfgOb, self.__configName)
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
            fastaPath = self.__cachePath = os.path.join(HERE, "test-output", "pharos-targets.fa")
            taxonPath = self.__cachePath = os.path.join(HERE, "test-output", "pharos-targets-taxon.tdd")
            ok = ptP.exportProteinFasta(fastaPath, taxonPath, addTaxonomy=True)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testPharosDataSite(self):
        try:
            url = "http://habanero.health.unm.edu/tcrd/download/latest.README"
            urlFallback = "https://unmtid-dbs.net/download/TCRD/latest.README"
            ok1 = ok2 = False
            try:
                res = requests.get(url, timeout=30)
                logger.info("Fetch url %r (response %r, len %r)", url, res.status_code, len(res.text))
                ok1 = (res.status_code == 200) and (len(res.text) > 500)
            except Exception:
                logger.error("Failed to fetch Pharos TCRD data from source %r", url)
            #
            res = requests.get(urlFallback, timeout=30)
            logger.info("Fetch fallback url %r (response %r, len %r)", urlFallback, res.status_code, len(res.text))
            ok2 = (res.status_code == 200) and (len(res.text) > 500)
            #
            self.assertTrue(ok1 or ok2)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def pharosTargetSuite():
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(PharosTargetProviderTests("testFetchPharosTargets"))
    # suiteSelect.addTest(PharosTargetProviderTests("testFetchAndLoadPharosTargets"))
    # suiteSelect.addTest(PharosTargetProviderTests("testExportPharosTargetFasta"))
    # suiteSelect.addTest(PharosTargetProviderTests("testExportPharosTargetFastaTax"))
    suiteSelect.addTest(PharosTargetProviderTests("testPharosDataSite"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = pharosTargetSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
