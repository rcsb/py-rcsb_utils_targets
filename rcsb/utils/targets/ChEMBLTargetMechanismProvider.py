##
#  File:           ChEMBLTargetMechanismProvider.py
#  Date:           9-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for ChEMBL target mechanism data.

"""

import datetime
import logging
import os.path
import time

from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

Settings.Instance().TIMEOUT = 10  # pylint: disable=no-member
Settings.Instance().MAX_LIMIT = 50  # pylint: disable=no-member
Settings.MAX_LIMIT = 50

logger = logging.getLogger(__name__)


class ChEMBLTargetMechanismProvider(StashableBase):
    """Accessors for ChEMBL target mechanism data."""

    def __init__(self, cachePath, useCache):
        #
        self.__cachePath = cachePath
        self.__dirName = "ChEMBL-target-mechanism"
        super(ChEMBLTargetMechanismProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        baseVersion = 28
        self.__version = baseVersion
        logger.info("ChEMBL API MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
        self.__aD = self.__reload(self.__dirPath, useCache)

    def testCache(self, minCount=0):
        if minCount == 0:
            return True
        if self.__aD and (len(self.__aD) > minCount):
            logger.info("Mechanism data for (%d) targets", len(self.__aD))
            return True
        return False

    def getAssignmentVersion(self):
        return self.__version

    def getTargetMechanismDataPath(self):
        return os.path.join(self.__dirPath, "chembl-target-mechanism-data.json")

    def __reload(self, dirPath, useCache):
        startTime = time.time()
        aD = {}
        fU = FileUtil()
        fU.mkdir(dirPath)
        targetMechanismFilePath = self.getTargetMechanismDataPath()
        #
        if useCache and fU.exists(targetMechanismFilePath):
            logger.info("useCache %r using %r", useCache, targetMechanismFilePath)
            qD = self.__mU.doImport(targetMechanismFilePath, fmt="json")
            aD = qD["mechanism"] if "mechanism" in qD else {}
        #
        logger.info("Completed reload of (%d) at %s (%.4f seconds)", len(aD), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        #
        return aD

    def getTargetMechanisms(self, targetChEMBLId):
        try:
            return self.__aD[targetChEMBLId] if targetChEMBLId in self.__aD else []
        except Exception:
            return []

    def hasTargetMechanism(self, targetChEMBLId):
        try:
            return targetChEMBLId in self.__aD
        except Exception:
            return False

    def fetchTargetMechanismData(self, targetChEMBLIdList, skipExisting=True, chunkSize=50):
        """Get cofactor mechanism data for the input ChEMBL target list.

        Args:
            targetChEMBLIdList (list): list of ChEMBL target identifiers
            skipExisting (bool, optional): reuse any existing cached data (default: True)
            chunkSize(int, optional): ChEMBL API batch size for fetches (default: 50)

        Returns:
          bool:  True for success or False otherwise

        """
        atL = [
            "action_type",
            "molecule_chembl_id",
            "action_type",
            "mechanism_of_action",
            "max_phase",
            "target_chembl_id",
        ]
        targetD = self.__aD if self.__aD else {}
        idList = []
        if skipExisting:
            for tId in targetChEMBLIdList:
                if tId in self.__aD:
                    continue
                idList.append(tId)
        else:
            idList = targetChEMBLIdList

        numToProcess = len(idList)
        logger.info("Fetching mechanism data for (%d/%d)", numToProcess, len(targetChEMBLIdList))
        ok = False
        try:
            for ii in range(0, len(idList), chunkSize):
                logger.info("Begin chunk at ii %d/%d", ii, numToProcess)
                mch = new_client.mechanism  # pylint: disable=no-member
                mch.set_format("json")
                mDL = mch.filter(target_chembl_id__in=idList[ii : ii + chunkSize]).only(atL)

                logger.info("Results (%d)", len(mDL))
                if mDL:
                    for mD in mDL:
                        targetD.setdefault(mD["target_chembl_id"], []).append(self.__mechanismSelect(atL, mD))
                #
                logger.info("Completed chunk starting at (%d)", ii)
                tS = datetime.datetime.now().isoformat()
                vS = datetime.datetime.now().strftime("%Y-%m-%d")
                ok = self.__mU.doExport(self.getTargetMechanismDataPath(), {"version": vS, "created": tS, "mechanism": targetD}, fmt="json", indent=3)
                logger.info("Wrote completed chunk starting at (%d) (%r)", ii, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __mechanismSelect(self, atL, aD):
        return {at: aD[at] if at in aD else None for at in atL}
