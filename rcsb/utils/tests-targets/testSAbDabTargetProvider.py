##
# File:    SAbDabTargetProviderTests.py
# Author:  J. Westbrook
# Date:    30-Nov-2020
# Version: 0.001
#
# Update:
#   19-Oct-2026 dwp  Add SAbDab2 assignment tests
#
#
##
"""
Tests for utilities managing SAbDab target data.
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

from rcsb.utils.targets.SAbDabTargetProvider import SAbDabTargetProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class SAbDabTargetProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__fastaPath = os.path.join(HERE, "test-output", "sabdab-targets.fa")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testFetchSAbDabTargets(self):
        stP = SAbDabTargetProvider(cachePath=self.__cachePath, useCache=False)
        ok = stP.testCache()
        self.assertTrue(ok)
        ok = stP.exportFasta(self.__fastaPath)
        self.assertTrue(ok)

    def testSAbDabAssignments(self):
        """Test SAbDab2 assignments (heavy/light chains, antigen names, light chain type, V gene subgroups)."""
        stP = SAbDabTargetProvider(cachePath=self.__cachePath, useCache=False)
        self.assertGreater(len(stP.getAssignments()), 20000)
        # 5BK0: light chain A has V gene subgroup (IMGT) and a Kappa chain type; lookup by short and extended ID
        for instId in ["5bk0.A", "pdb_00005bk0.A"]:
            self.assertTrue(stP.hasAssignment(instId))
            self.assertEqual("Kappa", stP.getAssignment(instId, "light_ctype"))
            self.assertEqual("IGKV2", stP.getAssignment(instId, "light_subclass"))
            self.assertIn("Circumsporozoite protein NANP 5-mer", stP.getAssignment(instId, "antigen_name"))
        # Heavy chain: subgroup, and no light chain-only keys
        self.assertEqual("IGHV4", stP.getAssignment("4jn2.H", "heavy_subclass"))
        self.assertIsNone(stP.getAssignment("4jn2.H", "light_ctype"))
        # Chains lacking IMGT V gene subgroups still get other assignments (e.g., 10FD)
        self.assertEqual("Kappa", stP.getAssignment("10fd.E", "light_ctype"))
        self.assertIsNone(stP.getAssignment("10fd.E", "light_subclass"))
        self.assertEqual("Glutamate receptor ionotropic, NMDA 2B", stP.getAssignment("10fd.E", "antigen_name"))
        # Parenthetical organisms are stripped
        for _, aD in stP.getAssignments().items():
            for ky in ["heavy_subclass", "light_subclass"]:
                self.assertNotIn("(", aD.get(ky, ""))
        # Ions and haptens in the source data are retained
        self.assertIn("CHLORIDE ION", stP.getAssignment("10bt.A", "antigen_name"))
        self.assertIn("GLYCEROL", stP.getAssignment("4jn2.H", "antigen_name"))
        self.assertFalse(stP.hasAssignment("0xxx.Z"))


def fetchSAbDabTargets():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SAbDabTargetProviderTests("testFetchSAbDabTargets"))
    suiteSelect.addTest(SAbDabTargetProviderTests("testSAbDabAssignments"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchSAbDabTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
