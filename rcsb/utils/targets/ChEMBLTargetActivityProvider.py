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
from chembl_webresource_client.settings import Settings

from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.io.UrlRequestUtil import UrlRequestUtil

Settings.Instance().TIMEOUT = 10  # pylint: disable=no-member
Settings.Instance().MAX_LIMIT = 50  # pylint: disable=no-member
Settings.MAX_LIMIT = 50

logger = logging.getLogger(__name__)


class ChEMBLTargetActivityWorker(object):
    """A skeleton worker class that implements the interface expected by the multiprocessing module
    for fetching ChEMBL activity data --
    """

    def __init__(self, **kwargs):
        _ = kwargs

    def fetchActivity(self, dataList, procName, optionsD, workingDir):
        """Fetch ChEMBL activity for the input ChEMBL target identifier list"""
        _ = workingDir
        successList = []
        failList = []
        retList = []
        diagList = []
        #
        try:
            chunkSize = optionsD.get("chunkSize", 50)
            atL = optionsD.get("attributeList", [])
            maxActivity = optionsD.get("maxActivity", None)

            for ii in range(0, len(dataList), chunkSize):
                logger.info("Begin chunk at ii %d/(for %d targets)", ii, len(dataList))
                try:
                    act = new_client.activity  # pylint: disable=no-member
                    act.set_format("json")
                    actDL = (
                        act.filter(target_chembl_id__in=dataList[ii : ii + chunkSize])
                        .filter(standard_type__in=["IC50", "Ki", "EC50", "Kd"])
                        .filter(standard_value__isnull=False)
                        .order_by("-standard_value")
                        .only(atL)
                    )
                    logger.info("Results (%d)", len(actDL))
                    if actDL:
                        for actD in actDL[:maxActivity] if maxActivity else actDL:
                            actD = self.__activitySelect(atL, actD)
                            actD["molecule_name"], actD["inchi_key"], _ = self.__getMoleculeDetails(actD["molecule_chembl_id"])
                            actD["action"], actD["moa"], actD["max_phase"] = self.__getMechanismDetails(actD["molecule_chembl_id"])
                            retList.append((actD["target_chembl_id"], actD))
                except Exception as e:
                    logger.exception("Failing for chunk starting at %d with %s", ii, str(e)[:200])
            successList = sorted(set(dataList) - set(failList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s built target interactions for %d/%d entries failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))
        #
        return successList, retList, diagList

    def __activitySelect(self, atL, aD):
        return {at: aD[at] if at in aD else None for at in atL}

    def __getMoleculeDetails(self, chemblId):
        name = inchiKey = smiles = None
        try:
            atL = ["pref_name", "molecule_structures"]
            molecule = new_client.molecule  # pylint: disable=no-member
            molecule.set_format("json")
            tD = molecule.filter(molecule_chembl_id__exact=chemblId).only(atL)[0]
            name = tD["pref_name"] if tD and "pref_name" in tD else None
            if tD and "molecule_structures" in tD and tD["molecule_structures"]:
                smiles = tD["molecule_structures"]["canonical_smiles"] if tD and "molecule_structures" in tD and "canonical_smiles" in tD["molecule_structures"] else None
                inchiKey = tD["molecule_structures"]["standard_inchi_key"] if tD and "molecule_structures" in tD and "standard_inchi_key" in tD["molecule_structures"] else None
        except Exception as e:
            logger.exception("Failing for %s with %s", chemblId, str(e))
        return name, inchiKey, smiles

    def __getMechanismDetails(self, chemblId):
        actionType = moa = maxPhase = None
        try:
            atL = ["action_type", "mechanism_of_action", "max_phase"]
            mechanism = new_client.mechanism  # pylint: disable=no-member
            mechanism.set_format("json")
            tD = mechanism.filter(molecule_chembl_id__exact=chemblId).only(atL)[0]
            actionType = tD["action_type"] if tD and "action_type" in tD else None
            moa = tD["mechanism_of_action"] if tD and "mechanism_of_action" in tD else None
            maxPhase = tD["max_phase"] if tD and "max_phase" in tD else None
        except Exception as e:
            logger.exception("Failing for %s with %s", chemblId, str(e))
        return actionType, moa, maxPhase


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
        self.__aD, self.__allIdD = self.__reload(self.__dirPath, useCache)

    def testCache(self, minCount=0):
        if minCount == 0:
            return True
        if self.__aD and (len(self.__aD) > minCount):
            logger.info("Activity data cached for (%d) targets", len(self.__aD))
            return True
        return False

    def getAssignmentVersion(self):
        return self.__version

    def getTargetActivityDataPath(self):
        return os.path.join(self.__dirPath, "chembl-target-activity-data.json")

    def reload(self):
        self.__aD, self.__allIdD = self.__reload(self.__dirPath, True)

    def __reload(self, dirPath, useCache):
        startTime = time.time()
        aD = {}
        allIdD = {}
        fU = FileUtil()
        fU.mkdir(dirPath)
        targetActivityFilePath = self.getTargetActivityDataPath()
        #
        if useCache and fU.exists(targetActivityFilePath):
            logger.info("useCache %r using %r", useCache, targetActivityFilePath)
            qD = self.__mU.doImport(targetActivityFilePath, fmt="json")
            aD = qD["activity"] if "activity" in qD else {}
            idL = qD["all_ids"] if "all_ids" in qD else []
            allIdD = {k: k in aD for k in idL}
        #
        logger.info(
            "Completed reload (%d activities) (%d tried identifiers) at %s (%.4f seconds)",
            len(aD),
            len(allIdD),
            time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
            time.time() - startTime,
        )
        #
        return aD, allIdD

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
                qCmtD = self.__decodeComment(queryId)
                tS = qCmtD["chemblId"]
                tL = tS.split(",")
                chemblIdList.extend(tL)
            chemblIdList = list(set(chemblIdList))
            logger.info("Total targets from sequence matching (%d)", len(chemblIdList))
        except Exception as e:
            logger.exception("Failing for %r with %s", sequenceMatchFilePath, str(e))
        return chemblIdList

    def fetchTargetActivityData(self, targetChEMBLIdList, skip="none", maxActivity=10, chunkSize=50):
        """Get cofactor activity data for the input ChEMBL target list.

        Args:
            targetChEMBLIdList (list): list of ChEMBL target identifiers
            skip (str, optional): Skip searching identifiers previously tried|matched|none (default: none)
            chunkSize (int, optional): ChEMBL API batch size for fetches (default: 50)
            maxActivity (int, optional): maximum number of activity records to return per target (default: 10)

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
        logger.info("Fetching activities for starting target list (%d) skip option %r", len(targetChEMBLIdList), skip)
        idList = []
        if skip == "matched":
            for tId in targetChEMBLIdList:
                if tId in self.__aD:
                    continue
                idList.append(tId)
        elif skip == "tried":
            for tId in targetChEMBLIdList:
                if tId in self.__allIdD:
                    continue
                idList.append(tId)
        else:
            idList = targetChEMBLIdList

        numToProcess = len(idList)
        logger.info("Filtered target list (%d)", len(idList))
        try:
            for ii in range(0, len(idList), chunkSize):
                logger.info("Begin chunk at ii %d/%d", ii, numToProcess)
                try:
                    act = new_client.activity  # pylint: disable=no-member
                    act.set_format("json")
                    actDL = (
                        act.filter(target_chembl_id__in=idList[ii : ii + chunkSize])
                        .filter(standard_type__in=["IC50", "Ki", "EC50", "Kd"])
                        .filter(standard_value__isnull=False)
                        .order_by("-standard_value")
                        .only(atL)
                    )
                    logger.info("Results for index %d (%d)", ii, len(actDL))
                    if actDL:
                        for actD in actDL[:maxActivity] if maxActivity else actDL:
                            actD = self.__activitySelect(atL, actD)
                            actD["molecule_name"], actD["inchi_key"], _ = self.getMoleculeDetails(actD["molecule_chembl_id"])
                            actD["action_type"], actD["moa"], actD["max_phase"] = self.getMechanismDetails(actD["molecule_chembl_id"])
                            targetD.setdefault(actD["target_chembl_id"], []).append(actD)
                    #
                except Exception as e:
                    logger.exception("Failing with chunk at index %d with %s", ii, str(e)[:200])

                logger.info("Completed chunk starting at (%d)", ii)
                tS = datetime.datetime.now().isoformat()
                vS = datetime.datetime.now().strftime("%Y-%m-%d")
                ok = self.__mU.doExport(self.getTargetActivityDataPath(), {"version": vS, "created": tS, "activity": targetD, "all_ids": idList}, fmt="json", indent=3)
                logger.info("Wrote completed chunk starting at (%d) (%r)", ii, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __activitySelect(self, atL, aD):
        return {at: aD[at] if at in aD else None for at in atL}

    def getMoleculeDetails(self, chemblId):
        name = inchiKey = smiles = None
        try:
            atL = ["pref_name", "molecule_structures"]
            molecule = new_client.molecule  # pylint: disable=no-member
            molecule.set_format("json")
            tD = molecule.filter(molecule_chembl_id__exact=chemblId).only(atL)[0]
            name = tD["pref_name"] if tD and "pref_name" in tD else None
            if tD and "molecule_structures" in tD and tD["molecule_structures"]:
                smiles = tD["molecule_structures"]["canonical_smiles"] if tD and "molecule_structures" in tD and "canonical_smiles" in tD["molecule_structures"] else None
                inchiKey = tD["molecule_structures"]["standard_inchi_key"] if tD and "molecule_structures" in tD and "standard_inchi_key" in tD["molecule_structures"] else None
        except Exception as e:
            logger.exception("Failing for %s with %s", chemblId, str(e))
        return name, inchiKey, smiles

    def getStatusDetails(self):
        try:
            baseUrl = "https://www.ebi.ac.uk"
            endPoint = "chembl/api/data/status.json"
            hL = []
            pD = {}
            ureq = UrlRequestUtil()
            ret, retCode = ureq.get(baseUrl, endPoint, pD, headers=hL, returnContentType="JSON")
            logger.info("retCode %r ret %r", retCode, ret)
            tS = ret["chembl_db_version"] if "chembl_db_version" in ret else None
            version = tS.split("_")[1] if tS and tS.split("_")[1] else None
            releaseDateString = ret["chembl_release_date"] if "chembl_release_date" in ret else None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return version, releaseDateString

    def getMechanismDetails(self, chemblId):
        actionType = moa = maxPhase = None
        try:
            atL = ["action_type", "mechanism_of_action", "max_phase"]
            mechanism = new_client.mechanism  # pylint: disable=no-member
            mechanism.set_format("json")
            tD = mechanism.filter(molecule_chembl_id__exact=chemblId).only(atL)[0]
            actionType = tD["action_type"] if tD and "action_type" in tD else None
            moa = tD["mechanism_of_action"] if tD and "mechanism_of_action" in tD else None
            maxPhase = tD["max_phase"] if tD and "max_phase" in tD else None
        except Exception as e:
            logger.exception("Failing for %s with %s", chemblId, str(e))
        return actionType, moa, maxPhase

    def fetchTargetActivityDataMulti(self, targetChEMBLIdList, skip="none", maxActivity=10, chunkSize=50, numProc=4):
        """Get cofactor activity data for the input ChEMBL target list (multiprocessing mode).

        Args:
            targetChEMBLIdList (list): list of ChEMBL target identifiers
            skip (str, optional): Skip searching identifiers previously tried|matched|none (default: none)
            chunkSize(int, optional): ChEMBL API outer batch size for fetches (default: 50)
            numProc (int, optional): number processes to invoke. (default: 4)
            maxActivity (int, optional): number of activity records to return per target. (default: 10)

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
        logger.info("Fetching activities for starting target list (%d) skip option %r", len(targetChEMBLIdList), skip)
        ok = False
        targetD = self.__aD if self.__aD else {}
        idList = []
        if skip == "matched":
            for tId in targetChEMBLIdList:
                if tId in self.__aD:
                    continue
                idList.append(tId)
        elif skip == "tried":
            for tId in targetChEMBLIdList:
                if tId in self.__allIdD:
                    continue
                idList.append(tId)
        else:
            idList = targetChEMBLIdList

        numToProcess = len(idList)
        allIdL = []
        logger.info("Filtered target list (%d)", numToProcess)
        ok = numToProcess == 0
        try:
            for ii in range(0, len(idList), chunkSize):
                logger.info("Begin outer chunk at ii %d (total targets %d)", ii, numToProcess)
                tIdList = idList[ii : ii + chunkSize]
                #
                tD = self.__getActivityMulti(tIdList, atL, maxActivity=maxActivity, numProc=numProc, chunkSize=5)
                targetD.update(tD)
                #
                allIdL.extend(tIdList)
                #
                logger.info("Completed outer chunk (%d) total processed targets (%d/%d)", ii, len(allIdL), numToProcess)
                tS = datetime.datetime.now().isoformat()
                vS = datetime.datetime.now().strftime("%Y-%m-%d")
                ok = self.__mU.doExport(self.getTargetActivityDataPath(), {"version": vS, "created": tS, "activity": targetD, "all_ids": allIdL}, fmt="json", indent=3)
                logger.info("Wrote completed outer chunk starting at (%d) (%r)", ii, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __getActivityMulti(self, idList, atL, maxActivity=None, numProc=2, chunkSize=5):
        """ """
        rD = {}
        ctaW = ChEMBLTargetActivityWorker()
        mpu = MultiProcUtil(verbose=True)
        optD = {"attributeList": atL, "chunkSize": chunkSize, "maxActivity": maxActivity}
        mpu.setOptions(optD)
        mpu.set(workerObj=ctaW, workerMethod="fetchActivity")
        ok, failList, resultList, _ = mpu.runMulti(dataList=idList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("Target Id activity failures (%d): %r", len(failList), failList)
        #
        for (targetId, actD) in resultList[0]:
            rD.setdefault(targetId, []).append(actD)
        #
        logger.info("Completed with multi-proc status %r failures %r total targets with data (%d)", ok, len(failList), len(rD))
        return rD

    def __decodeComment(self, comment, separator="|"):
        dD = {}
        try:
            ti = iter(comment.split(separator))
            dD = {tup[1]: tup[0] for tup in zip(ti, ti)}
        except Exception:
            pass
        return dD
