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
from difflib import SequenceMatcher
import logging
import os.path
import time

from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from yaml import SequenceEndEvent

Settings.Instance().TIMEOUT = 10  # pylint: disable=no-member
Settings.Instance().MAX_LIMIT = 50  # pylint: disable=no-member
Settings.MAX_LIMIT = 50

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
        logger.info("ChEMBL API MAX_LIMIT %r", Settings.Instance().MAX_LIMIT)  # pylint: disable=no-member
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

    def getTargetIdList(self, sequenceMatchFilePath):
        chemblIdList = []
        try:
            mD = self.__mU.doImport(sequenceMatchFilePath, fmt="json")
            # --- cofactor list
            chemblIdList = []
            for queryId in mD:
                tS = queryId.split("|")[2]
                tL = tS.split(",")
                chemblIdList.extend(tL)
            chemblIdList = list(set(chemblIdList))
            logger.info("Total matched targets (%d)", len(chemblIdList))
        except Exception as e:
            logger.exception("Failing for %r with %s", sequenceMatchFilePath, str(e))
        return chemblIdList

    def fetchTargetActivityData(self, targetChEMBLIdList, skipExisting=True, chunkSize=50):
        """Get cofactor activity data for the input ChEMBL target list.

        Args:
            targetChEMBLIdList (list): list of ChEMBL target identifiers
            skipExisting (bool, optional): reuse any existing cached data (default: True)
            chunkSize(int, optional): ChEMBL API batch size for fetches (default: 50)

        Returns:
          bool:  True for success or False otherwise

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
        ok = False
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
        return ok

    def __activitySelect(self, atL, aD):

        return {at: aD[at] if at in aD else None for at in atL}
