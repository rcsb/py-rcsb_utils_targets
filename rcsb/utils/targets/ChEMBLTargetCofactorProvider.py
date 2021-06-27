##
#  File:           ChEMBLTargetCofactorProvider.py
#  Date:           15-Jun-2021 jdw
#
#  Updated:
#
##
"""
Accessors for ChEMBL target cofactors.
"""

import datetime
import logging
import math
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.targets.ChEMBLTargetActivityProvider import ChEMBLTargetActivityProvider
from rcsb.utils.targets.ChEMBLTargetProvider import ChEMBLTargetProvider

logger = logging.getLogger(__name__)


class ChEMBLTargetCofactorProvider(StashableBase):
    """Accessors for ChEMBL target cofactors."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirName = "ChEMBL-cofactors"
        super(ChEMBLTargetCofactorProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fD = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=1):
        logger.info("ChEMBL cofactor count %d", len(self.__fD["cofactors"]) if "cofactors" in self.__fD else 0)
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
        return os.path.join(self.__dirPath, "ChEMBL-cofactor-data.json")

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

            Example activity record -

                        "CHEMBL3243": [
                    {
                        "assay_chembl_id": "CHEMBL655768",
                        "assay_description": "In vitro inhibitory activity against recombinant human CD45 using fluorescein diphosphate (FDP) as a substrate",
                        "assay_type": "B",
                        "canonical_smiles": "COC(=O)c1ccc(C2=CC(=O)C(=O)c3ccccc32)cc1",
                        "ligand_efficiency": {
                        "bei": "19.78",
                        "le": "0.36",
                        "lle": "3.11",
                        "sei": "9.57"
                        },
                        "molecule_chembl_id": "CHEMBL301254",
                        "parent_molecule_chembl_id": "CHEMBL301254",
                        "pchembl_value": "5.78",
                        "standard_relation": "=",
                        "standard_type": "IC50",
                        "standard_units": "nM",
                        "standard_value": "1650.0",
                        "target_chembl_id": "CHEMBL3243"
                    },
        """
        rDL = []
        mD = self.__mU.doImport(sequenceMatchFilePath, fmt="json")
        #
        chP = ChEMBLTargetProvider(cachePath=self.__cachePath, useCache=False)
        # --- cofactor list
        chemblIdList = []
        for queryId, matchDL in mD.items():
            tS = queryId.split("|")[2]
            tL = tS.split(",")
            chemblIdList.extend(tL)
        chemblIdList = list(set(chemblIdList))
        logger.info("Total cofactors for matched targets (%d)", len(chemblIdList))
        # ---
        #
        chaP = ChEMBLTargetActivityProvider(cachePath=self.__cachePath, useCache=True)
        #
        provenanceSource = "ChEMBL"
        refScheme = "PDB entity"
        assignVersion = chP.getAssignmentVersion()
        for queryId, matchDL in mD.items():
            unpId = queryId.split("|")[0].strip()
            queryTaxId = queryId.split("|")[4]
            chemblIdL = queryId.split("|")[2].split(",")
            if queryTaxId == "-1":
                logger.info("Skipping target %r (%r)", unpId, chemblIdL)
                continue
            queryName = chP.getTargetDescription(unpId)
            for chemblId in chemblIdL:
                if not chaP.hasTargetActivity(chemblId):
                    logger.debug("Skipping target %r (%r)", unpId, chemblId)
                    continue
                #
                for matchD in matchDL:
                    tL = matchD["target"].split("|")
                    entryId = tL[0].split("_")[0]
                    entityId = tL[0].split("_")[1]
                    #
                    taDL = chaP.getTargetActivity(chemblId)
                    logger.debug("Target %r has (%d) activity records", chemblId, len(taDL))
                    # ---
                    actL = []
                    for taD in taDL:
                        if taD["assay_type"] in ["B", "F"]:
                            if taD["standard_units"] == "nM":
                                try:
                                    pV = -math.log10(float(taD["standard_value"]) * 10.0e-9)
                                    actD = {
                                        "cofactor_id": taD["molecule_chembl_id"],
                                        "assay_id": taD["assay_chembl_id"],
                                        "assay_description": taD["assay_description"],
                                        "measurement_type": "p" + taD["standard_type"],
                                        "measurement_value": round(pV, 2),
                                        "smiles": taD["canonical_smiles"],
                                        "molecule_name": taD["molecule_name"],
                                        "inchi_key": taD["inchi_key"],
                                        "action": taD["action"],
                                        "moa": taD["moa"],
                                        "max_phase": taD["max_phase"],
                                    }

                                    actL.append(actD)
                                except Exception as e:
                                    logger.exception("Failing for tAD %r with %s", taD, str(e))

                    # ---
                    actL = self.__activityListSelect(actL, crmpObj, maxActivity=maxActivity)
                    # ---
                    rD = {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "query_uniprot_id": unpId,
                        "query_id": chemblId,
                        "query_id_type": "ChEMBL",
                        "query_name": queryName,
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
        # vS = datetime.datetime.now().strftime("%Y-%m-%d")
        vS = assignVersion
        ok = self.__mU.doExport(fp, {"version": vS, "created": tS, "cofactors": qD}, fmt="json", indent=3)
        return ok

    def __activityListSelect(self, activityDL, crmpOb=None, maxActivity=5):
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
                            activityD["prd_id"] = localId
                        else:
                            activityD["chem_comp_id"] = localId
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
