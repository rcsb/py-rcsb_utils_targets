##
#  File:           ChEMBLTargetActivityProvider.py
#  Date:           9-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for ChEMBL target activity data.

"""

import datetime
import logging
import os.path
import time

from chembl_webresource_client.new_client import new_client

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase


logger = logging.getLogger(__name__)


class ChEMBLTargetActivityProvider(StashableBase):
    """Accessors for ChEMBL target activity data."""

    def __init__(self, cachePath, useCache):
        #
        self.__cachePath = cachePath
        self.__dirName = "ChEMBL-target-activity"
        super(ChEMBLTargetActivityProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        baseVersion = 28
        self.__version = baseVersion
        self.__aD = self.__reload(self.__dirPath, useCache)

    def testCache(self, minCount=0):
        if minCount == 0:
            return True
        if self.__aD and (len(self.__aD) > minCount):
            logger.info("Activity data for (%d) targets", len(self.__aD))
            return True
        return False

    def getAssignmentVersion(self):
        return self.__version

    def getTargetActivityDataPath(self):
        return os.path.join(self.__dirPath, "chembl-target-activity-data.json")

    def __reload(self, dirPath, useCache):
        startTime = time.time()
        aD = {}
        fU = FileUtil()
        fU.mkdir(dirPath)
        targetActivityFilePath = self.getTargetActivityDataPath()
        #
        if useCache and fU.exists(targetActivityFilePath):
            logger.info("useCache %r using %r", useCache, targetActivityFilePath)
            qD = self.__mU.doImport(targetActivityFilePath, fmt="json")
            aD = qD["activity"] if "activity" in qD else {}
        #
        logger.info("Completed reload of (%d) at %s (%.4f seconds)", len(aD), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        #
        return aD

    def getTargetActivity(self, targetChEMBLId):
        try:
            return self.__aD[targetChEMBLId] if targetChEMBLId in self.__aD else []
        except Exception:
            return []

    def hasTargetActivity(self, targetChEMBLId):
        try:
            return targetChEMBLId in self.__aD
        except Exception:
            return False

    def fetchActivityData(self, targetChEMBLIdList, skipExisting=True, chunkSize=50):
        """Get cofactor activity data for the input ChEMBL target list.

        Args:
            targetChEMBLIdList (list): list of ChEMBL target identifiers

        Returns:
          (dict, dict):  {targetChEMBId: {activity data}}, {moleculeChEMBId: {activity data}}

        """
        atL = [
            "assay_chembl_id",
            "assay_description",
            "assay_type",
            "canonical_smiles",
            "ligand_efficiency",
            "molecule_chembl_id",
            "parent_molecule_chembl_id",
            "pchembl_value",
            "standard_relation",
            "standard_type",
            "standard_units",
            "standard_value",
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
        try:
            for ii in range(0, len(idList), chunkSize):
                logger.info("Begin chunk at ii %d/%d", ii, numToProcess)
                act = new_client.activity  # pylint: disable=no-member
                act.set_format("json")
                actDL = (
                    act.filter(target_chembl_id__in=idList[ii : ii + chunkSize]).filter(standard_type__in=["IC50", "Ki", "EC50", "Kd"]).filter(standard_value__isnull=False).only(atL)
                )
                logger.info("Results (%d)", len(actDL))
                if actDL:
                    for actD in actDL:
                        targetD.setdefault(actD["target_chembl_id"], []).append(self.__activitySelect(atL, actD))
                #
                logger.info("Completed chunk starting at (%d)", ii)
                tS = datetime.datetime.now().isoformat()
                vS = datetime.datetime.now().strftime("%Y-%m-%d")
                ok = self.__mU.doExport(self.getTargetActivityDataPath(), {"version": vS, "created": tS, "activity": targetD}, fmt="json", indent=3)
                logger.info("Wrote completed chunk starting at (%d) (%r)", ii, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return targetD

    def __activitySelect(self, atL, aD):

        return {at: aD[at] if at in aD else None for at in atL}

    def getMechanismData(self, targetChEMBLIdList):
        """Get mechanism data for the input ChEMBL target list.

        Args:
            targetChEMBLIdList (list): list of ChEMBL target identifiers

        Returns:
          (dict):  dictionary  {ChEMBId: {mechanism data}}

        """
        oD = {}
        chunkSize = 50
        try:
            for ii in range(0, len(targetChEMBLIdList), chunkSize):
                mch = new_client.mechanism  # pylint: disable=no-member
                mch.set_format("json")
                mDL = mch.filter(target_chembl_id__in=targetChEMBLIdList[ii : ii + chunkSize])
                if mDL:
                    logger.info("mDL (%d)", len(mDL))
                    for mD in mDL:
                        oD.setdefault(mD["target_chembl_id"], []).append(mD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD
