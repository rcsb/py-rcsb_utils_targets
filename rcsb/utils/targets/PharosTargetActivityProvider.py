##
#  File:           PharosTargetActivityProvider.py
#  Date:           17-Jun-2021 jdw
#
#  Updated:
#
##
"""
Accessors for Pharos target activity data.

"""

import datetime
import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase


logger = logging.getLogger(__name__)


class PharosTargetActivityProvider(StashableBase):
    """Accessors for Pharos target activity data."""

    def __init__(self, cachePath, useCache):
        #
        self.__cachePath = cachePath
        self.__dirName = "Pharos-target-activity"
        super(PharosTargetActivityProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        #
        self.__aD, self.__tD, self.__version = self.__reload(self.__dirPath, useCache)

    def testCache(self, minCount=0):
        if minCount == 0:
            return True
        if self.__aD and (len(self.__aD) > minCount):
            logger.info("Cached Pharos activity data for (%d) targets", len(self.__aD))
            return True
        return False

    def getAssignmentVersion(self):
        return self.__version

    def getTargetActivityDataPath(self):
        return os.path.join(self.__dirPath, "pharos-target-activity-data.json")

    def __reload(self, dirPath, useCache):
        startTime = time.time()
        aD = {}
        tD = {}
        version = None
        fU = FileUtil()
        fU.mkdir(dirPath)
        targetActivityFilePath = self.getTargetActivityDataPath()
        #
        if useCache and fU.exists(targetActivityFilePath):
            logger.info("useCache %r using %r", useCache, targetActivityFilePath)
            qD = self.__mU.doImport(targetActivityFilePath, fmt="json")
            aD = qD["activity"] if "activity" in qD else {}
            tD = qD["targets"] if "targets" in qD else {}
            version = qD["version"] if "version" in qD else None
        #
        logger.info("Completed reload of (%d) activity records at %s (%.4f seconds)", len(aD), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        #
        return aD, tD, version

    def getTargetActivity(self, pharosTargetId):
        try:
            return self.__aD[pharosTargetId] if pharosTargetId in self.__aD else []
        except Exception:
            return []

    def hasTargetActivity(self, pharosTargetId):
        try:
            return pharosTargetId in self.__aD
        except Exception:
            return False

    def hasTargetInfo(self, pharosTargetId):
        try:
            return pharosTargetId in self.__tD
        except Exception:
            return False

    def getTargetInfo(self, pharosTargetId, ky):
        try:
            return self.__tD[pharosTargetId][ky]
        except Exception:
            return None

    def fetchTargetActivityData(self):
        targetD = {}
        cofactorFilePath = os.path.join(self.__cachePath, "Pharos-targets", "drug_activity.tdd")
        cfDL = self.__mU.doImport(cofactorFilePath, fmt="tdd", rowFormat="dict")
        targetD = self.__extactCofactorData(cfDL)
        #
        cofactorFilePath = os.path.join(self.__cachePath, "Pharos-targets", "cmpd_activity.tdd")
        cfDL = self.__mU.doImport(cofactorFilePath, fmt="tdd", rowFormat="dict")
        targetD.update(self.__extactCofactorData(cfDL))
        #
        targetFilePath = os.path.join(self.__cachePath, "Pharos-targets", "protein.tdd")
        tDL = self.__mU.doImport(targetFilePath, fmt="tdd", rowFormat="dict")
        targetDetailsD = self.__getTargetDetails(tDL)
        #
        pharosReadmePath = os.path.join(self.__cachePath, "Pharos-targets", "pharos-readme.txt")
        readmeLines = self.__mU.doImport(pharosReadmePath, fmt="list")
        self.__version = readmeLines[0].split(" ")[1][1:] if readmeLines else "6"
        #
        tS = datetime.datetime.now().isoformat()
        # vS = datetime.datetime.now().strftime("%Y-%m-%d")
        vS = self.__version
        ok = self.__mU.doExport(self.getTargetActivityDataPath(), {"version": vS, "created": tS, "activity": targetD, "targets": targetDetailsD}, fmt="json", indent=3)
        return ok

    def __extactCofactorData(self, cfDL):
        """Extract ids, activity and moa data for drugs and cofactors from the Pharos schema dump files.

        Args:
            cfDL (list): list of dictionaries of containing pharos exported db data.

        Returns:
            dict: dictionary of extracted cofactor data
        """
        try:
            qD = {}
            targetD = {}
            dupD = {}
            for cfD in cfDL:
                tId = cfD["target_id"]
                qD = {}
                qD["smiles"] = cfD["smiles"] if "smiles" in cfD and cfD["smiles"] not in ["N", "NULL"] else None
                qD["chemblId"] = cfD["cmpd_chemblid"] if "cmpd_chemblid" in cfD else None
                qD["chemblId"] = cfD["cmpd_id_in_src"] if "catype" in cfD and cfD["catype"].upper() == "CHEMBL" else qD["chemblId"]
                qD["chemblId"] = qD["chemblId"] if qD["chemblId"] not in ["N", "NULL"] else None
                qD["pubChemId"] = cfD["cmpd_pubchem_cid"] if "cmpd_pubchem_cid" in cfD and cfD["cmpd_pubchem_cid"] not in ["NULL"] else None
                #
                qD["activity"] = cfD["act_value"] if "act_value" in cfD and cfD["act_value"] != "NULL" else None
                qD["activityType"] = cfD["act_type"] if "act_type" in cfD and cfD["act_type"] != "NULL" else None
                if qD["activity"] is not None:
                    qD["activity"] = float(qD["activity"])
                    # qD["activityUnits"] = "nM"
                #
                qD["action"] = cfD["action"] if "action" in cfD and cfD["moa"] == "1" else None
                qD["pharmacology"] = cfD["nlm_drug_info"] if "nlm_drug_info" in cfD else None
                tS = cfD["drug"] if "drug" in cfD else None
                tS = cfD["cmpd_name_in_src"] if "cmpd_name_in_src" in cfD and cfD["cmpd_name_in_src"] != "NULL" else tS
                #
                if tS and tS.startswith("US"):
                    tSL = tS.split(",")
                    qD["patents"] = [t for t in tSL if t.startswith("US")]
                elif tS:
                    # handle the double colon :: separators
                    tSL = tS.split("::")
                    if not tSL:
                        qD["molecule_name"] = tS
                    else:
                        for tS in tSL:
                            if tS.startswith("CHEMBL"):
                                continue
                            elif tS.startswith("US"):
                                ttL = tS.split(",")
                                qD.setdefault("patents", []).append(ttL[0])
                            else:
                                qD["molecule_name"] = tS

                #
                pmId = None
                tS = cfD["reference"] if "reference" in cfD else None
                if tS and "pubmed" in tS:
                    pmId = tS.split("/")[-1]
                tS = cfD["pubmed_ids"].split(",")[0] if "pubmed_ids" in cfD and cfD["pubmed_ids"] else pmId
                qD["pubmedId"] = tS if tS and tS not in ["NULL"] else None
                #
                if qD["activity"] and qD["chemblId"] and (qD["chemblId"], qD["activityType"], qD["action"]) not in dupD:
                    dupD[(qD["chemblId"], qD["activityType"], qD["action"])] = True
                    targetD.setdefault(tId, []).append({ky: qD[ky] for ky in qD if qD[ky] is not None})

            #
        except Exception as e:
            logger.exception("Failing with %r %s", qD, str(e))
        return targetD

    def __getTargetDetails(self, targetDL):
        # Pharos protein target - protein.tdd
        # id	name	description	uniprot	up_version	geneid	sym	family	chr	seq	dtoid	stringid	dtoclass
        rD = {}
        try:
            for targetD in targetDL:
                proteinId = targetD["id"]
                unpId = targetD["uniprot"] if "uniprot" in targetD and targetD["uniprot"] != "NULL" else None
                descr = targetD["description"] if "description" in targetD and targetD["description"] != "NULL" else None
                geneId = targetD["geneid"] if "geneid" in targetD and targetD["geneid"] != "NULL" else None
                dtoId = targetD["dtoid"] if "dtoid" in targetD and targetD["dtoid"] != "NULL" else None
                dtoClass = targetD["dtoclass"] if "dtoclass" in targetD and targetD["dtoclass"] != "NULL" else None
                rD[proteinId] = {"unpId": unpId, "name": descr, "geneId": geneId, "dtoId": dtoId, "dtoClass": dtoClass}
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return rD
