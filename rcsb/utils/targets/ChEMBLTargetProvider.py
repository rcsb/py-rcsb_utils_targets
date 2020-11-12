##
#  File:           ChEMBLTargetProvider.py
#  Date:           9-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for ChEMBL target assignments.

"""

import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ChEMBLTargetProvider:
    """Accessors for ChEMBL target assignments."""

    def __init__(self, **kwargs):
        #
        self.__version = "0.50"
        cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirPath = os.path.join(cachePath, "CACHE", "ChEMBL")

        chemblDbUrl = kwargs.get("ChEMBLDbUrl", "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__phD = self.__reload(chemblDbUrl, self.__dirPath, useCache=useCache)
        #

    def testCache(self):
        return True

    def __reload(self, chemblDbUrl, dirPath, useCache=True):
        startTime = time.time()
        ok = False
        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        baseVersion = 27
        # ChEMBL current version 27,...
        # template:  chembl_27.fa.gz
        #
        targetFileName = "chembl_" + str(baseVersion) + ".fa.gz"
        mappingFileName = "chembl_uniprot_mapping.txt"
        #
        chemblTargetPath = os.path.join(dirPath, targetFileName)
        chemblMappingPath = os.path.join(dirPath, mappingFileName)
        #
        if useCache and self.__mU.exists(chemblMappingPath):
            logger.info("useCache %r using %r and %r", useCache, chemblTargetPath, chemblMappingPath)
            ok = True
        else:
            url = os.path.join(chemblDbUrl, mappingFileName)
            ok = fU.get(url, chemblMappingPath)
            logger.info("Fetching url %s path %s", url, chemblMappingPath)
            #
            for vers in range(baseVersion, baseVersion + 10):
                targetFileName = "chembl_" + str(vers) + ".fa.gz"
                chemblTargetPath = os.path.join(dirPath, targetFileName)
                url = os.path.join(chemblDbUrl, targetFileName)
                ok = fU.get(url, chemblTargetPath)
                logger.info("Fetching url %s path %s", url, chemblTargetPath)
                if ok:
                    break
            #
            logger.info("Completed fetches at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            #

        return ok
