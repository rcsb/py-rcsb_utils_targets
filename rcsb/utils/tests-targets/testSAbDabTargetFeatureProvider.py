##
# File:    SAbDabTargetFeatureProviderTests.py
# Author:  J. Westbrook
# Date:    11-Jun-2021
# Version: 0.001
#
# Update:
#   20-Oct-2026 dwp  Look up assignments by extended ID and entity features without forcing case;
#                    check the feature data version is the assignment data version
#
#
##
"""
Tests for utilities managing SAbDab target feature data.
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

from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.targets.SAbDabTargetFeatureProvider import SAbDabTargetFeatureProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class SAbDabTargetFeatureProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dataPath = os.path.join(HERE, "test-data")
        #
        self.__seqMatchResultsPath = os.path.join(HERE, "test-data", "sabdab-vs-pdbprent-filtered-results.json.gz")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildSAbDabTargetsFeatures(self):
        stfP = SAbDabTargetFeatureProvider(cachePath=self.__cachePath, useCache=False)
        ok = stfP.testCache()
        self.assertFalse(ok)
        ok = stfP.buildFeatureList(self.__seqMatchResultsPath)
        self.assertTrue(ok)
        stfP = SAbDabTargetFeatureProvider(cachePath=self.__cachePath, useCache=True)
        ok = stfP.testCache()
        self.assertTrue(ok)
        ok = stfP.hasFeatures("4G6K_2")
        self.assertTrue(ok)
        fL = stfP.getFeatures("4G6K_2")
        self.assertGreaterEqual(len(fL), 3)
        # Assignments are keyed by extended ID
        ok = stfP.hasAssignment("pdb_00005bk0.A")
        self.assertTrue(ok)
        ok = stfP.getAssignment("pdb_00005bk0.A", "antigen_name")
        self.assertTrue(ok)
        # Feature data version is the SAbDab assignment data version (also used for the entity features)
        self.assertRegex(stfP.getVersion(), r"^\d{4}-\d{2}-\d{2}$")
        self.assertEqual(stfP.getVersion(), fL[0]["assignment_version"])

    @unittest.skip("Bootstrap test")
    def testSabdabCacheBootstrap(self):
        try:
            stfP = SAbDabTargetFeatureProvider(cachePath=self.__cachePath, useCache=True)
            ok = stfP.testCache()
            self.assertTrue(ok)
            configPath = os.path.join(self.__dataPath, "sabdab-config.yml")
            configName = "site_info_remote_configuration"
            cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
            ok = stfP.backup(cfgOb, configName, useGit=True, useStash=True)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def buildSAbDabFeaturesTargets():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SAbDabTargetFeatureProviderTests("testBuildSAbDabTargetsFeatures"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = buildSAbDabFeaturesTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
