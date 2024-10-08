##
#  File:           DrugBankTargetCofactorProvider.py
#  Date:           12-Jun-2021 jdw
#
#  Updated:
#  20-Aug-2024 dwp Add support for loading and accessing data on MongoDB
##
"""
Accessors for DrugBank target cofactors.
"""

import datetime
import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.targets.DrugBankTargetProvider import DrugBankTargetProvider
from rcsb.utils.targets.TargetCofactorDbProvider import TargetCofactorDbProvider

logger = logging.getLogger(__name__)


class DrugBankTargetCofactorProvider(StashableBase):
    """Accessors for DrugBank target cofactors."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__useCache = kwargs.get("useCache", True)
        self.__fmt = kwargs.get("fmt", "pickle")
        self.__dirName = "DrugBank-cofactors"
        self.__cofactorResourceName = "drugbank"
        super(DrugBankTargetCofactorProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fD = self.__reload(self.__dirPath, self.__useCache, self.__fmt)
        #

    def testCache(self, minCount=590):
        logger.info("DrugBank cofactor count %d", len(self.__fD["cofactors"]) if "cofactors" in self.__fD else 0)
        if self.__fD and "cofactors" in self.__fD and len(self.__fD["cofactors"]) > minCount:
            return True
        else:
            return False

    def hasTarget(self, rcsbEntityId):
        return rcsbEntityId.upper() in self.__fD["cofactors"]

    def getTargets(self, rcsbEntityId):
        try:
            return self.__fD["cofactors"][rcsbEntityId.upper()]
        except Exception:
            return []

    def __getCofactorDataPath(self, fmt="json"):
        fExt = "json" if fmt == "json" else "pic"
        return os.path.join(self.__dirPath, "drugbank-cofactor-data.%s" % fExt)

    def reload(self):
        self.__fD = self.__reload(self.__dirPath, useCache=True, fmt=self.__fmt)
        return True

    def __reload(self, dirPath, useCache, fmt):
        startTime = time.time()
        fD = {}

        ok = False
        cofactorPath = self.__getCofactorDataPath(fmt=fmt)
        #
        logger.info("useCache %r cofactorPath %r", useCache, cofactorPath)
        if useCache and self.__mU.exists(cofactorPath):
            fD = self.__mU.doImport(cofactorPath, fmt=fmt)
        else:
            fU = FileUtil()
            fU.mkdir(dirPath)
        # ---
        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return fD

    def getCofactorDataDict(self):
        return self.__fD["cofactors"]

    def buildCofactorList(self, sequenceMatchFilePath, crmpObj=None, lnmpObj=None):
        """Build target cofactor list for the matching entities in the input sequence match file.

        Args:
            sequenceMatchFilePath (str): sequence match output file path
            crmpObj (obj, optional): instance of ChemRefMappingProviderObj(). Defaults to None
            lnmpObj (obj, optional): instance of LigandNeighborMappingProviderObj(). Defaults to None.

        Returns:
            bool: True for success or False otherwise
        """
        rDL = []
        dbP = DrugBankTargetProvider(cachePath=self.__cachePath, useCache=True)
        if not dbP.testCache():
            logger.warning("Skipping build of target cofactor list because DrugBank Target data is missing.")
            return False
        mD = self.__mU.doImport(sequenceMatchFilePath, fmt="json")
        #
        provenanceSource = "DrugBank"
        refScheme = "PDB entity"
        assignVersion = str(dbP.getAssignmentVersion())
        for queryId, matchDL in mD.items():
            qCmtD = self.__decodeComment(queryId)
            unpId = qCmtD["uniprotId"]
            queryTaxId = qCmtD["taxId"] if "taxId" in qCmtD else None
            if not dbP.hasCofactor(unpId) or queryTaxId == "-1":
                logger.info("Skipping target %r", unpId)
                continue
            #
            # --
            chemCompNeighborsD = {}
            if lnmpObj:
                for matchD in matchDL:
                    tCmtD = self.__decodeComment(matchD["target"])
                    entryId = tCmtD["entityId"].split("_")[0]
                    entityId = tCmtD["entityId"].split("_")[1]
                    rcsbEntityId = entryId + "_" + entityId
                    chemCompIdList = lnmpObj.getLigandNeighbors(rcsbEntityId)
                    chemCompNeighborsD.update({k: True for k in chemCompIdList})
            # --
            #
            for matchD in matchDL:
                tCmtD = self.__decodeComment(matchD["target"])
                entryId = tCmtD["entityId"].split("_")[0]
                entityId = tCmtD["entityId"].split("_")[1]
                # --
                dbDL = dbP.getCofactors(unpId)
                # --
                cfDL = []
                for dbD in dbDL:
                    cfD = {}
                    cfD["cofactor_id"] = dbD["drugbank_id"]
                    cfD["molecule_name"] = dbD["name"]
                    cfD["target_name"] = dbD["target_name"]
                    # cfD["description"] = dbD["description"]
                    cfD["moa"] = dbD["moa"]
                    # cfD["pharmacology"] = dbD["pharmacology"]
                    cfD["inchi_key"] = dbD["inchi_key"]
                    cfD["smiles"] = dbD["smiles"]
                    cfD["pubmed_ids"] = dbD["pubmed_ids"]
                    cfD = self.__addLocalIds(cfD, crmpObj)
                    #
                    if "chem_comp_id" in cfD and cfD["chem_comp_id"] in chemCompNeighborsD:
                        cfD["neighbor_in_pdb"] = "Y"
                    else:
                        cfD["neighbor_in_pdb"] = "N"
                    #
                    cfDL.append(cfD)
                # ---
                queryName = cfDL[0]["target_name"] if cfDL and "target_name" in cfDL[0] else None
                # ---
                # aligned_target.entity_beg_seq_id (current target is PDB entity in json)
                # aligned_target.target_beg_seq_id (current query is target seq in json)
                # aligned_target.length
                fpL = []
                if "alignedRegions" in matchD:
                    fpL = [
                        {
                            "entity_beg_seq_id": arD["targetBegin"],
                            "target_beg_seq_id": arD["queryBegin"],
                            "length": arD["targetEnd"] - arD["targetBegin"],
                        }
                        for arD in matchD["alignedRegions"]
                    ]
                else:
                    fpL = [
                        {
                            "entity_beg_seq_id": matchD["targetBegin"],
                            "target_beg_seq_id": matchD["queryBegin"],
                            "length": matchD["alignLen"],
                        }
                    ]
                # ---
                rD = {
                    "entry_id": entryId,
                    "entity_id": entityId,
                    "query_uniprot_id": unpId,
                    "query_id": unpId,
                    "query_id_type": "DrugBank",
                    "query_name": queryName,
                    "provenance_source": provenanceSource,
                    "reference_scheme": refScheme,
                    "assignment_version": assignVersion,
                    "query_taxonomy_id": int(queryTaxId) if queryTaxId else None,
                    "target_taxonomy_id": int(matchD["targetTaxId"]) if "targetTaxId" in matchD else None,
                    "aligned_target": fpL,
                    "taxonomy_match_status": matchD["taxonomyMatchStatus"] if "taxonomyMatchStatus" in matchD else None,
                    "lca_taxonomy_id": matchD["lcaTaxId"] if "lcaTaxId" in matchD else None,
                    "lca_taxonomy_name": matchD["lcaTaxName"] if "lcaTaxName" in matchD else None,
                    "lca_taxonomy_rank": matchD["lcaRank"] if "lcaRank" in matchD else None,
                    "cofactors": cfDL,
                }
                rDL.append(rD)
        #
        qD = {}
        for rD in rDL:
            eId = rD["entry_id"] + "_" + rD["entity_id"]
            qD.setdefault(eId, []).append(rD)
        fp = self.__getCofactorDataPath(fmt=self.__fmt)
        tS = datetime.datetime.now().isoformat()
        # vS = datetime.datetime.now().strftime("%Y-%m-%d")
        vS = assignVersion
        ok = self.__mU.doExport(fp, {"version": vS, "created": tS, "cofactors": qD}, fmt=self.__fmt, indent=3)
        return ok

    def __addLocalIds(self, cfD, crmpOb=None):
        #
        if crmpOb:
            localIdL = crmpOb.getLocalIds("DRUGBANK", cfD["cofactor_id"])
            if localIdL:
                localId = localIdL[0]
                if localId.startswith("PRD_"):
                    cfD["prd_id"] = localId
                else:
                    cfD["chem_comp_id"] = localId
        return cfD

    def __decodeComment(self, comment, separator="|"):
        dD = {}
        try:
            ti = iter(comment.split(separator))
            dD = {tup[1]: tup[0] for tup in zip(ti, ti)}
        except Exception:
            pass
        return dD

    def loadCofactorData(self, cfgOb, **kwargs):
        """Load cofactor data to MongoDB.
        """
        ok = False
        try:
            tcDbP = TargetCofactorDbProvider(
                cachePath=self.__cachePath,
                cfgOb=cfgOb,
                cofactorResourceName=self.__cofactorResourceName,
                **kwargs
            )
            ok = tcDbP.loadCofactorData(self.__cofactorResourceName, self)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        logger.info("%r cofactor data DB load status (%r)", self.__cofactorResourceName, ok)
        return ok


class DrugBankTargetCofactorAccessor:
    def __init__(self, cachePath, cfgOb=None, **kwargs):
        """
        Accessor class for fetching cofactor data from MongoDB.
        """
        self.__cofactorResourceName = "drugbank"
        self.__cfgOb = cfgOb
        self.__cachePath = cachePath
        #
        self.__tcDbP = TargetCofactorDbProvider(
            cachePath=self.__cachePath,
            cfgOb=self.__cfgOb,
            cofactorResourceName=self.__cofactorResourceName,
            **kwargs
        )

    def testCache(self, minCount=0):
        docCount = self.__tcDbP.cofactorDbCount()
        logger.info("Loaded %s cofactor DB count %d", self.__cofactorResourceName, docCount)
        return docCount >= minCount

    def reload(self):
        return True

    def getTargets(self, rcsbEntityId, dataFieldName="rcsb_cofactors", **kwargs):
        return self.__tcDbP.fetchCofactorData(rcsbEntityId, dataFieldName, **kwargs)
