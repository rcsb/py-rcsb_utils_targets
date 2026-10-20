##
# File:    SAbDabTargetProviderTests.py
# Author:  J. Westbrook
# Date:    30-Nov-2020
# Version: 0.001
#
# Update:
#   19-Oct-2026 dwp  Add SAbDab2 assignment tests
#   20-Oct-2026 dwp  Add primary source availability and fallback tests (incl. gzipped backup and version date)
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
import gzip
import json
import shutil
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

    def testPrimarySourcesAvailable(self):
        """Fail if either primary data source is unavailable or invalid (fallbacks are allowed in production, but not here)."""
        stP = SAbDabTargetProvider(cachePath=os.path.join(HERE, "test-output", "CACHE-primary"), useCache=False, allowFallback=False)
        sD = stP.getSources()
        self.assertFalse(sD["thera"]["fallback"])
        self.assertFalse(sD["assignments"]["fallback"])
        self.assertEqual(time.strftime("%Y-%m-%d"), stP.getAssignmentVersion())
        self.assertTrue(stP.testCache())
        self.assertGreater(len(stP.getAssignments()), 20000)

    def testStashBackupSources(self):
        """Test the default stash backup sources (fetched over the network) by making the primary sources unavailable."""
        missing = os.path.join(HERE, "test-output", "does-not-exist")
        stP = SAbDabTargetProvider(
            cachePath=os.path.join(HERE, "test-output", "CACHE-stash"), useCache=False, targetUrl=missing + ".csv", assignmentUrl=missing + ".json"
        )
        sD = stP.getSources()
        self.assertTrue(sD["thera"]["fallback"])
        self.assertTrue(sD["assignments"]["fallback"])
        self.assertTrue(sD["assignments"]["url"].endswith("sabdab_summary_all.json.gz"))
        self.assertRegex(stP.getAssignmentVersion(), r"^\d{4}-\d{2}-\d{2}$")
        self.assertTrue(stP.testCache())
        self.assertGreater(len(stP.getAssignments()), 20000)
        self.assertEqual("Kappa", stP.getAssignment("5bk0.A", "light_ctype"))

    def testFallbackSources(self):
        """Test fallback to backup sources (local files) when the primary sources are missing or invalid."""
        dirPath = os.path.join(HERE, "test-output", "fallback")
        shutil.rmtree(dirPath, ignore_errors=True)
        os.makedirs(dirPath)
        theraPath = os.path.join(dirPath, "thera-backup.csv")
        with open(theraPath, "w", encoding="utf-8") as ofh:
            ofh.write("Therapeutic,Format,CH1 Isotype,VD LC,Highest_Clin_Trial (Feb '25),Est. Status,Target,Conditions Approved,Conditions Active,HeavySequence,LightSequence\n")
            ofh.write("Fooumab,Whole mAb,IgG1,Kappa,Phase-III,Active,FOO,na,Cancer,EVQLVESGGG,DIQMTQSPSS\n")
        assignPath = os.path.join(dirPath, "assign-backup.json.gz")
        rD = {"PDB ID": "pdb_00001abc", "model": 0, "antigen_instances": [{"antigen name": "FOO", "entity type": "PROTEIN"}], "light chain": {"PDB auth_asym_id": "L", "chain type": "κ", "V gene subgroup": None}}
        with gzip.open(assignPath, "wt", encoding="utf-8") as ofh:
            json.dump({"version": "2026-01-02", "records": [rD]}, ofh, ensure_ascii=False)
        # An HTML page served with a 200 status is an invalid primary source
        htmlPath = os.path.join(dirPath, "page.html")
        with open(htmlPath, "w", encoding="utf-8") as ofh:
            ofh.write("<!doctype html><html></html>")
        kw = {"cachePath": os.path.join(dirPath, "CACHE"), "useCache": False, "minAssignmentRecords": 1, "targetFallbackUrl": theraPath, "assignmentFallbackUrl": assignPath}
        for primaryTheraUrl, primaryAssignUrl in [(os.path.join(dirPath, "missing.csv"), os.path.join(dirPath, "missing.json")), (htmlPath, htmlPath)]:
            stP = SAbDabTargetProvider(targetUrl=primaryTheraUrl, assignmentUrl=primaryAssignUrl, **kw)
            sD = stP.getSources()
            self.assertTrue(sD["thera"]["fallback"])
            self.assertTrue(sD["assignments"]["fallback"])
            self.assertEqual("Kappa", stP.getAssignment("1abc.L", "light_ctype"))
            self.assertEqual("FOO", stP.getAssignment("pdb_00001abc.L", "antigen_name"))
            # Version is the date recorded in the backup file, not today's date
            self.assertEqual("2026-01-02", stP.getAssignmentVersion())
            self.assertEqual("2026-01-02", stP.getSources()["assignments"]["version"])
        # Fallbacks can be disabled
        with self.assertRaises(ValueError):
            SAbDabTargetProvider(targetUrl=htmlPath, assignmentUrl=htmlPath, allowFallback=False, **{k: v for k, v in kw.items() if k != "useCache"}, useCache=False)
        # Too few records in a source fails validation (and with no usable fallback, raises)
        with self.assertRaises(ValueError):
            SAbDabTargetProvider(targetUrl=theraPath, assignmentUrl=assignPath, targetFallbackUrl=theraPath, assignmentFallbackUrl=assignPath, cachePath=os.path.join(dirPath, "CACHE2"), useCache=False)

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
    suiteSelect.addTest(SAbDabTargetProviderTests("testPrimarySourcesAvailable"))
    suiteSelect.addTest(SAbDabTargetProviderTests("testStashBackupSources"))
    suiteSelect.addTest(SAbDabTargetProviderTests("testFallbackSources"))
    suiteSelect.addTest(SAbDabTargetProviderTests("testSAbDabAssignments"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchSAbDabTargets()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
