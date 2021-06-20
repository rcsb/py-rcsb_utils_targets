##
#  File:           PharosTargetCofactorProvider.py
#  Date:           15-Jun-2021 jdw
#
#  Updated:
#
##
"""
Accessors for Pharos target cofactor data.
"""

import datetime
import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.targets.PharosTargetActivityProvider import PharosTargetActivityProvider

logger = logging.getLogger(__name__)


class PharosTargetCofactorProvider(StashableBase):
    """Accessors for Pharos target cofactors."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirName = "Pharos-cofactors"
        super(PharosTargetCofactorProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fD = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=1):
        logger.info("Pharos cofactor count %d", len(self.__fD["cofactors"]) if "cofactors" in self.__fD else 0)
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
        return os.path.join(self.__dirPath, "Pharos-cofactor-data.json")

    def __reload(self, dirPath, **kwargs):
        startTime = time.time()
        fD = {}
        useCache = kwargs.get("useCache", True)
        ok = False
        cofactorPath = self.__getCofactorDataPath()
        #
        logger.info("useCache %r cofactorPath %r", useCache, cofactorPath)
        if useCache and self.__mU.exists(cofactorPath):
            fD = self.__mU.doImport(cofactorPath, fmt="json")
            ok = True
        else:
            fU = FileUtil()
            fU.mkdir(dirPath)
        # ---
        logger.info("Completed reload with status (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return fD

    def buildCofactorList(self, sequenceMatchFilePath, crmpObj=None, maxActivity=5):
        """Build target cofactor list for the matching entities in the input sequence match file.

        Args:
            sequenceMatchFilePath (str): sequence match output file path
            crmpObj (obj, optional): instance of ChemRefMappingProviderObj()
            maxActivity (int, optional): maximum number of prioritized activity records per target

        Returns:
            bool: True for success or False otherwise

            Example Pharos activity record -

            {
            "version": "2021-06-17",
            "created": "2021-06-17T11:10:54.563394",
            "activity": {
                "2232": [
                    {
                        "smiles": "CC(=CCC\\\\C(=C/Cc1c(O)cc(O)c(C(=O)CCc2ccc(O)cc2)c1O)\\\\C)C",
                        "chemblId": "CHEMBL3360923",
                        "pubChemId": "118724585",
                        "activity": 6.0,
                        "activityType": "IC50",
                        "activityUnits": "nM",
                        "name": "1-[3-(3,7-dimethylocta-2,6-dien-1-yl)-2,4,6-trihydroxyphenyl]-3-(4-hydroxyphenyl)propan-1-one",
                        "pubmedId": "25375026",
                        "patent": "USxxxxxx",
                    }, ...
        """
        rDL = []
        mD = self.__mU.doImport(sequenceMatchFilePath, fmt="json")
        #
        # --- cofactor list
        chemblIdList = []
        for queryId, matchDL in mD.items():
            tS = queryId.split("|")[2]
            tL = tS.split(",")
            chemblIdList.extend(tL)
        chemblIdList = list(set(chemblIdList))
        logger.info("Total cofactors for matched targets (%d)", len(chemblIdList))
        # ---
        chaP = PharosTargetActivityProvider(cachePath=self.__cachePath, useCache=True)
        #
        provenanceSource = "Pharos"
        refScheme = "PDB entity"
        assignVersion = chaP.getAssignmentVersion()
        for queryId, matchDL in mD.items():
            # "O43508|uniprotId|7987|proteinId|9606|taxId"
            unpId = queryId.split("|")[0]
            queryTaxId = queryId.split("|")[4]
            pharosId = queryId.split("|")[2]
            if queryTaxId == "-1":
                logger.debug("Skipping target %r (%r)", unpId, pharosId)
                continue
            #
            if not chaP.hasTargetActivity(pharosId):
                logger.debug("Skipping target %r (%r)", unpId, pharosId)
                continue
            #
            for matchD in matchDL:
                tL = matchD["target"].split("|")
                entryId = tL[0].split("_")[0]
                entityId = tL[0].split("_")[1]
                #
                taDL = chaP.getTargetActivity(pharosId)
                logger.debug("Target %r has (%d) activity records", pharosId, len(taDL))
                actL = []
                cfDL = []
                chD = {}
                for taD in taDL:
                    if taD["chemblId"] in chD:
                        chD[taD["chemblId"]] = True
                        continue
                    # JDW leaving this in place but this is too large to include
                    cfD = {"cofactor_id": taD["chemblId"]}
                    cfDL.append(self.__addLocalIds(cfD, crmpObj))

                    actD = {
                        "cofactor_id": taD["chemblId"],
                        "cofactor_name": taD["name"] if "name" in taD else None,
                        "assay_description": None,
                        "measurement_type": "p" + taD["activityType"],
                        "measurement_value": taD["activity"],
                        "pubmed_id": taD["pubmedId"] if "pubmedId" in taD else None,
                        "patent_no": taD["patent"] if "patent" in taD else None,
                    }
                    actL.append(actD)
                #
                actL = self.__activityListSelect(actL, crmpObj, maxActivity=maxActivity)
                # ---
                rD = {
                    "entry_id": entryId,
                    "entity_id": entityId,
                    "query_uniprot_id": unpId,
                    "query_id": pharosId,
                    "query_id_type": "Pharos",
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
                    "cofactors": actL,
                }
                rDL.append(rD)
        #
        qD = {}
        for rD in rDL:
            eId = rD["entry_id"] + "_" + rD["entity_id"]
            qD.setdefault(eId, []).append(rD)
        #
        fp = self.__getCofactorDataPath()
        tS = datetime.datetime.now().isoformat()
        vS = datetime.datetime.now().strftime("%Y-%m-%d")
        ok = self.__mU.doExport(fp, {"version": vS, "created": tS, "cofactors": qD}, fmt="json", indent=3)
        return ok

    def __addLocalIds(self, cfD, crmpOb=None):
        #
        if crmpOb:
            localIdL = crmpOb.getLocalIds("CHEMBL", cfD["cofactor_id"])
            if localIdL:
                localId = localIdL[0]
                if localId.startswith("PRD_"):
                    cfD["prd_id"] = localId
                else:
                    cfD["chem_comp_id"] = localId
        return cfD

    def __activityListSelect(self, activityDL, crmpOb=None, maxActivity=5):
        """Prioritizing the activity data for locally mapped ligands and best binding examples.

        Args:
            activityDL (list): full list of activity objects
            crmpOb (obj, optional): instance of ChemRefMappingProvider(). Defaults to None.
            maxCount (int, optional): maximum number of activity object returned. Defaults to 5.

        Returns:
            list: prioritized and trimmed list of activity objects
        """
        retL = []
        mappedL = []
        unmappedL = activityDL
        if crmpOb:
            unmappedL = []
            # Select out the any cases for molecules that map to a local CC or BIRD
            for activityD in activityDL:
                localIdL = crmpOb.getLocalIds("CHEMBL", activityD["cofactor_id"])
                if localIdL:
                    # add the local identifiers
                    for localId in localIdL:
                        if localId.startswith("PRD_"):
                            activityD.setdefault("prd_id", []).append(localId)
                        else:
                            activityD.setdefault("chem_comp_id", []).append(localId)
                    mappedL.append(activityD)
                else:
                    unmappedL.append(activityD)
        #
        numLeft = maxActivity - len(mappedL)
        if numLeft > 0:
            unmappedL = sorted(unmappedL, key=lambda k: k["measurement_value"], reverse=True)
            retL = mappedL
            retL.extend(unmappedL[:numLeft])
            retL = sorted(retL, key=lambda k: k["measurement_value"], reverse=True)
        return retL
