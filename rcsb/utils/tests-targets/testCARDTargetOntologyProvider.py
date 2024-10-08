##
# File:    testCARDTargetOntologyProvider.py
# Author:  D. Piehl
# Date:    14-Mar-2023
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing CARD ontology data.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.targets.CARDTargetOntologyProvider import CARDTargetOntologyProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class CARDTargetOntologyProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testFetchCARDTargetOntology(self):
        coP = CARDTargetOntologyProvider(cachePath=self.__cachePath, useCache=False)
        ok = coP.buildOntologyData()
        self.assertTrue(ok)
        ok = coP.reload()
        self.assertTrue(ok)
        ok = coP.testCache()
        self.assertTrue(ok)
        #
        coP = None
        coP = CARDTargetOntologyProvider(cachePath=self.__cachePath, useCache=True)
        ok = coP.testCache()
        self.assertTrue(ok)
        aroIdL = ["ARO:3000096", "ARO:3000015", "ARO:3001110", "ARO:0000041", "ARO:0000039", "ARO:3000454"]
        for aroId in aroIdL:
            lineageL = coP.getLineage(aroId)
            logger.info("Lineage list for aroId %s: %r", aroId, lineageL)
            self.assertGreaterEqual(len(lineageL), 2)
        tnL = coP.getTreeNodeList()
        self.assertGreater(len(tnL), 100)
        aroIdL += ["ARO:1000003", "ARO:3000000", "ARO:0000076", "ARO:1000002", "ARO:3000708", "ARO:3000082", "ARO:3000045"]
        for treeNode in tnL:
            if treeNode["id"] in aroIdL:
                logger.info("tree node aroId %s: %r", aroId, treeNode)


def fetchCARDTargetOntology():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CARDTargetOntologyProviderTests("testFetchCARDTargetOntology"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchCARDTargetOntology()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
