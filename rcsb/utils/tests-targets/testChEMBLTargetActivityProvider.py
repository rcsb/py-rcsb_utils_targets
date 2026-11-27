##
# File:    ChEMBLTargetActivityProviderTests.py
# Author:  J. Westbrook
# Date:    15-Jun-2021
#
# Update:
#  27-Nov-2026 dwp  Add offline tests: targets in a failed ChEMBL API chunk are reported as failures and not
#                   recorded as tried
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
import tempfile
import unittest
from unittest import mock

from rcsb.utils.targets.ChEMBLTargetActivityProvider import ChEMBLTargetActivityProvider, ChEMBLTargetActivityWorker
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

    def testFetchStatus(self):
        try:
            ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=False)
            version, releaseDateString = ctP.getStatusDetails()
            if not version:
                logger.warning("status service request failed")
            else:
                self.assertGreaterEqual(int(version), 28)
                self.assertTrue(releaseDateString is not None)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchMoleculeData(self):
        try:
            ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=False)
            chemblId = "CHEMBL1421"
            name, inchiKey, _ = ctP.getMoleculeDetails(chemblId)
            self.assertEqual(name, "DASATINIB ANHYDROUS")
            self.assertEqual(inchiKey, "ZBNZXTGUTAYRHI-UHFFFAOYSA-N")
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchActivityDataMulti(self):
        try:
            ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=False)
            ok = ctP.testCache()
            self.assertFalse(ok)
            #
            tL = ["CHEMBL1987", "CHEMBL3243"]
            ok = ctP.fetchTargetActivityDataMulti(tL, numProc=2)
            self.assertTrue(ok)
            ctP.reload()
            ok = ctP.testCache()
            self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testFetchActivityData(self):
        try:
            ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=False)
            ok = ctP.testCache()
            self.assertFalse(ok)

            mL = ["CHEMBL200117"]
            for chemblId in mL:
                molD = ctP.getMoleculeDetails(chemblId)
                logger.info("%s molD %r", chemblId, molD)
            #
            tL = ["CHEMBL1987", "CHEMBL3243"]
            ok = ctP.fetchTargetActivityData(tL)
            self.assertTrue(ok)
            ctP.reload()
            ok = ctP.testCache()
            self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


class _FakeActivityQuery(object):
    """Stand-in for a chembl_webresource_client activity query: 'CHEMBL_BAD' fails on the first request,
    'CHEMBL_PARTIAL' fails after returning one record, and any other target returns one record."""

    def __init__(self):
        self.__idL = []

    def set_format(self, _fmt):
        pass

    def filter(self, **kwargs):
        if "target_chembl_id__in" in kwargs:
            self.__idL = list(kwargs["target_chembl_id__in"])
        return self

    def order_by(self, *_args):
        return self

    def only(self, _atL):
        return self

    def __records(self):
        for tId in self.__idL:
            if tId == "CHEMBL_PARTIAL":
                yield {"target_chembl_id": tId, "molecule_chembl_id": "CHEMBL1"}
                raise RuntimeError("Simulated failure while paging results")
            yield {"target_chembl_id": tId, "molecule_chembl_id": "CHEMBL1"}

    def __len__(self):
        if "CHEMBL_BAD" in self.__idL:
            raise RuntimeError("Simulated HTTP 500 from the ChEMBL API")
        return len(self.__idL)

    def __getitem__(self, sl):
        return self.__records()

    def __iter__(self):
        return self.__records()


class ChEMBLTargetActivityFailureTests(unittest.TestCase):
    """Offline tests (no ChEMBL API calls) for the handling of failed activity fetches."""

    def setUp(self):
        self.__tmpDir = tempfile.TemporaryDirectory()
        self.__cachePath = self.__tmpDir.name
        self.__mU = MarshalUtil(workPath=self.__cachePath)

    def tearDown(self):
        self.__tmpDir.cleanup()

    def __fakeClient(self):
        fakeClient = mock.MagicMock()
        fakeClient.activity = _FakeActivityQuery()
        return fakeClient

    def testWorkerReportsFailedChunks(self):
        """A failed chunk's targets are reported as failures and its partial results are dropped."""
        optD = {"attributeList": ["target_chembl_id", "molecule_chembl_id"], "chunkSize": 1, "maxActivity": 10}
        dataList = ["CHEMBL_OK", "CHEMBL_BAD", "CHEMBL_PARTIAL"]
        with mock.patch("rcsb.utils.targets.ChEMBLTargetActivityProvider.new_client", self.__fakeClient(), create=True):
            successList, retList, _ = ChEMBLTargetActivityWorker().fetchActivity(dataList, "test-proc", optD, self.__cachePath)
        self.assertEqual(successList, ["CHEMBL_OK"])
        self.assertEqual(sorted({tId for tId, _ in retList}), ["CHEMBL_OK"])

    def testFailedTargetsNotRecordedAsTried(self):
        """Failed targets are left out of 'all_ids', so the next 'tried' run fetches them again."""
        actD = {"target_chembl_id": "CHEMBL_OK", "molecule_chembl_id": "CHEMBL1"}
        ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=False)
        with mock.patch.object(ChEMBLTargetActivityProvider, "_ChEMBLTargetActivityProvider__getActivityMulti", return_value=({"CHEMBL_OK": [actD]}, ["CHEMBL_BAD"])):
            ok = ctP.fetchTargetActivityDataMulti(["CHEMBL_OK", "CHEMBL_BAD"], skip="tried", numProc=1)
        self.assertTrue(ok)
        qD = self.__mU.doImport(ctP.getTargetActivityDataPath(), fmt="json")
        self.assertEqual(sorted(qD["all_ids"]), ["CHEMBL_OK"])
        self.assertEqual(sorted(qD["activity"].keys()), ["CHEMBL_OK"])
        #
        # The next run tries only the failed target
        ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=True)
        retryActD = {"target_chembl_id": "CHEMBL_BAD", "molecule_chembl_id": "CHEMBL2"}
        with mock.patch.object(ChEMBLTargetActivityProvider, "_ChEMBLTargetActivityProvider__getActivityMulti", return_value=({"CHEMBL_BAD": [retryActD]}, [])) as mockMulti:
            ok = ctP.fetchTargetActivityDataMulti(["CHEMBL_OK", "CHEMBL_BAD"], skip="tried", numProc=1)
        self.assertTrue(ok)
        self.assertEqual(mockMulti.call_args[0][0], ["CHEMBL_BAD"])
        qD = self.__mU.doImport(ctP.getTargetActivityDataPath(), fmt="json")
        self.assertEqual(sorted(qD["all_ids"]), ["CHEMBL_BAD", "CHEMBL_OK"])
        self.assertEqual(sorted(qD["activity"].keys()), ["CHEMBL_BAD", "CHEMBL_OK"])

    def testSingleProcessRecordsOnlyCompletedChunks(self):
        """fetchTargetActivityData() records only the targets of chunks that completed."""
        ctP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=False)
        with mock.patch("rcsb.utils.targets.ChEMBLTargetActivityProvider.new_client", self.__fakeClient(), create=True):
            ok = ctP.fetchTargetActivityData(["CHEMBL_OK", "CHEMBL_BAD"], chunkSize=1)
        self.assertTrue(ok)
        qD = self.__mU.doImport(ctP.getTargetActivityDataPath(), fmt="json")
        self.assertEqual(sorted(qD["all_ids"]), ["CHEMBL_OK"])


def targetActivityFailureSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChEMBLTargetActivityFailureTests("testWorkerReportsFailedChunks"))
    suiteSelect.addTest(ChEMBLTargetActivityFailureTests("testFailedTargetsNotRecordedAsTried"))
    suiteSelect.addTest(ChEMBLTargetActivityFailureTests("testSingleProcessRecordsOnlyCompletedChunks"))
    return suiteSelect


def targetActivitySuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChEMBLTargetActivityProviderTests("testFetchStatus"))
    suiteSelect.addTest(ChEMBLTargetActivityProviderTests("testFetchMoleculeData"))
    suiteSelect.addTest(ChEMBLTargetActivityProviderTests("testFetchActivityData"))
    suiteSelect.addTest(ChEMBLTargetActivityProviderTests("testFetchActivityDataMulti"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = targetActivitySuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    mySuite = targetActivityFailureSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
