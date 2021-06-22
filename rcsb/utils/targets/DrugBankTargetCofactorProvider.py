##
#  File:           DrugBankTargetCofactorProvider.py
#  Date:           12-Jun-2021 jdw
#
#  Updated:
#
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

logger = logging.getLogger(__name__)


class DrugBankTargetCofactorProvider(StashableBase):
    """Accessors for DrugBank target cofactors."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirName = "Drugbank-cofactors"
        super(DrugBankTargetCofactorProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fD = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=590):
        logger.info("Drugbank feature count %d", len(self.__fD["cofactors"]) if "cofactors" in self.__fD else 0)
        if self.__fD and "cofactors" in self.__fD and len(self.__fD["cofactors"]) > minCount:
            return True
        else:
            return False

    def hasCofactor(self, rcsbEntityId):
        return rcsbEntityId.upper() in self.__fD["cofactors"]

    def getCofactors(self, rcsbEntityId):
        try:
            return self.__fD["cofactors"][rcsbEntityId.upper()]
        except Exception:
            return None

    def __getCofactorDataPath(self):
        return os.path.join(self.__dirPath, "DrugBank-cofactor-data.json")

    def __reload(self, dirPath, **kwargs):
        startTime = time.time()
        fD = {}
        useCache = kwargs.get("useCache", True)
        ok = False
        cofactorPath = self.__getCofactorDataPath()
        #
        logger.info("useCache %r featurePath %r", useCache, cofactorPath)
        if useCache and self.__mU.exists(cofactorPath):
            fD = self.__mU.doImport(cofactorPath, fmt="json")
        else:
            fU = FileUtil()
            fU.mkdir(dirPath)
        # ---
        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return fD

    def buildCofactorList(
        self,
        sequenceMatchFilePath,
        crmpObj=None,
    ):
        """Build target cofactor list for the matching entities in the input sequence match file.

        Args:
            sequenceMatchFilePath (str): sequence match output file path
            crmpObj (obj, optional): instance of ChemRefMappingProviderObj()

        Returns:
            bool: True for success or False otherwise
        """
        rDL = []
        dbP = DrugBankTargetProvider(cachePath=self.__cachePath, useCache=False)
        mD = self.__mU.doImport(sequenceMatchFilePath, fmt="json")
        #
        provenanceSource = "Drugbank"
        refScheme = "PDB entity"
        assignVersion = dbP.getAssignmentVersion()
        for queryId, matchDL in mD.items():
            unpId = queryId.split("|")[0]
            queryTaxId = queryId.split("|")[2].strip()
            if not dbP.hasCofactor(unpId) or queryTaxId == "-1":
                logger.info("Skipping target %r", unpId)
                continue
            for matchD in matchDL:
                tL = matchD["target"].split("|")
                entryId = tL[0].split("_")[0]
                entityId = tL[0].split("_")[1]
                # --
                dbDL = dbP.getCofactors(unpId)
                # --
                cfDL = []
                for dbD in dbDL:
                    cfD = {}
                    cfD["cofactor_id"] = dbD["drugbank_id"]
                    cfD["molecule_name"] = dbD["name"]
                    cfD["description"] = dbD["description"]
                    cfD["moa"] = dbD["moa"]
                    cfD["pharmacology"] = dbD["pharmacology"]
                    cfD["inchi_key"] = dbD["inchi_key"]
                    cfD["smiles"] = dbD["smiles"]
                    cfD["pubmed_ids"] = dbD["pubmed_ids"]
                    cfDL.append(self.__addLocalIds(cfD, crmpObj))
                # --
                rD = {
                    "entry_id": entryId,
                    "entity_id": entityId,
                    "query_uniprot_id": unpId,
                    "query_id": unpId,
                    "query_id_type": "DrugBank",
                    "provenance_source": provenanceSource,
                    "reference_scheme": refScheme,
                    "assignment_version": assignVersion,
                    "query_taxonomy_id": int(queryTaxId),
                    "target_taxonomy_id": int(matchD["targetTaxId"]),
                    "target_beg_seq_id": matchD["targetStart"],
                    "query_beg_seq_id": matchD["queryStart"],
                    "align_length": matchD["alignLen"],
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
        fp = self.__getCofactorDataPath()
        tS = datetime.datetime.now().isoformat()
        vS = datetime.datetime.now().strftime("%Y-%m-%d")
        ok = self.__mU.doExport(fp, {"version": vS, "created": tS, "cofactors": qD}, fmt="json", indent=3)
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
