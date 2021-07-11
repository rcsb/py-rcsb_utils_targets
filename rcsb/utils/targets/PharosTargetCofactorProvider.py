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
        logger.info("Pharos cached cofactor count %d", len(self.__fD["cofactors"]) if "cofactors" in self.__fD else 0)
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

    def reload(self):
        self.__fD = self.__reload(self.__dirPath, useCache=True)
        return True

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
        numCofactors = len(fD["cofactors"]) if fD and "cofactors" in fD else 0
        logger.info(
            "Completed reload of (%d) cofactors with status (%r) at %s (%.4f seconds)", numCofactors, ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime
        )
        return fD

    def buildCofactorList(self, sequenceMatchFilePath, crmpObj=None, lnmpObj=None, maxActivity=5):
        """Build target cofactor list for the matching entities in the input sequence match file.

        Args:
            sequenceMatchFilePath (str): sequence match output file path
            crmpObj (obj, optional): instance of ChemRefMappingProviderObj(). Defaults to None.
            lnmpObj (obj, optional): instance of LigandNeighborMappingProviderObj(). Defaults to None.
            maxActivity (int, optional): maximum number of prioritized activity records per target. Defaults to 5.

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
        # ---
        chaP = PharosTargetActivityProvider(cachePath=self.__cachePath, useCache=True)
        #
        provenanceSource = "Pharos"
        refScheme = "PDB entity"
        assignVersion = chaP.getAssignmentVersion()
        for queryId, matchDL in mD.items():
            # "O43508|uniprotId|7987|proteinId|9606|taxId"
            qCmtD = self.__decodeComment(queryId)
            unpId = qCmtD["uniprotId"]
            queryTaxId = qCmtD["taxId"] if "taxId" in qCmtD else None
            pharosId = qCmtD["proteinId"]
            if queryTaxId == "-1":
                logger.debug("Skipping target with missing taxonomy %r (%r)", unpId, pharosId)
                continue
            #
            if not chaP.hasTargetActivity(pharosId):
                logger.debug("Skipping target with no activities %r (%r)", unpId, pharosId)
                # continue
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
            queryName = chaP.getTargetInfo(pharosId, "name")
            # --
            for matchD in matchDL:
                tCmtD = self.__decodeComment(matchD["target"])
                entryId = tCmtD["entityId"].split("_")[0]
                entityId = tCmtD["entityId"].split("_")[1]
                rcsbEntityId = entryId + "_" + entityId
                #
                taDL = chaP.getTargetActivity(pharosId)
                logger.debug("Target %r has (%d) activity records", pharosId, len(taDL))
                actL = []
                # cfDL = []
                chD = {}
                for taD in taDL:
                    if taD["chemblId"] in chD:
                        chD[taD["chemblId"]] = True
                        continue

                    actD = {
                        "cofactor_id": taD["chemblId"],
                        "cofactor_name": taD["molecule_name"] if "name" in taD else None,
                        "measurement_type": "p" + taD["activityType"],
                        "measurement_value": taD["activity"],
                        "pubmed_ids": [taD["pubmedId"]] if "pubmedId" in taD else None,
                        "patent_nos": taD["patents"] if "patents" in taD else None,
                        "smiles": taD["smiles"] if "smiles" in taD else None,
                        "action": taD["action"] if "action" in taD else None,
                        "pharmacology": taD["pharmacology"] if "pharmacology" in taD else None,
                    }
                    actD = self.__addLocalIds(actD, crmpObj=crmpObj)
                    actL.append(actD)
                #
                actL = self.__activityListSelect(actL, chemCompNeighborsD, maxActivity=maxActivity)
                if not actL:
                    logger.debug("No Pharos cofactors for %s %s", pharosId, unpId)
                # ---
                if "alignedRegions" in matchD:
                    tBegSeqId = ",".join([str(arD["targetBegin"]) for arD in matchD["alignedRegions"]])
                    qBegSeqId = ",".join([str(arD["queryBegin"]) for arD in matchD["alignedRegions"]])
                    alignLen = ",".join([str(arD["targetEnd"] - arD["targetBegin"]) for arD in matchD["alignedRegions"]])
                else:
                    tBegSeqId = matchD["targetStart"]
                    qBegSeqId = matchD["queryStart"]
                    alignLen = matchD["alignLen"]
                # ---
                rD = {
                    "entry_id": entryId,
                    "entity_id": entityId,
                    "query_uniprot_id": unpId,
                    "query_id": pharosId,
                    "query_id_type": "Pharos",
                    "query_name": queryName,
                    "provenance_source": provenanceSource,
                    "reference_scheme": refScheme,
                    "assignment_version": assignVersion,
                    "query_taxonomy_id": int(queryTaxId) if queryTaxId else None,
                    "target_taxonomy_id": int(matchD["targetTaxId"]) if "targetTaxId" in matchD else None,
                    #
                    "target_beg_seq_id": tBegSeqId,
                    "query_beg_seq_id": qBegSeqId,
                    "align_length": alignLen,
                    #
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

    def __addLocalIds(self, cfD, crmpObj=None):
        #
        if crmpObj:
            localIdL = crmpObj.getLocalIds("CHEMBL", cfD["cofactor_id"])
            if localIdL:
                localId = localIdL[0]
                if localId.startswith("PRD_"):
                    cfD["prd_id"] = localId
                else:
                    cfD["chem_comp_id"] = localId
        return cfD

    def __activityListSelect(self, activityDL, chemCompNeighborsD, maxActivity=5):
        """Prioritizing the activity data for locally mapped neighbor ligands and the best binding examples.

        Args:
            activityDL (list): full list of activity objects
            chemCompNeighborsD (dict, optional): index of all chemical components with neighbor interactions to the query target. Defaults {}.
            maxCount (int, optional): maximum number of activity object returned. Defaults to 5.

        Returns:
            list: prioritized and trimmed list of activity objects
        """
        retL = []
        mappedNeighborL = []
        unmappedL = activityDL

        if chemCompNeighborsD:
            unmappedL = []
            # Select out the any cases for molecules that map to a neighbor chemical component.
            for activityD in activityDL:
                if "chem_comp_id" in activityD and activityD["chem_comp_id"] in chemCompNeighborsD:
                    activityD["neighbor_in_pdb"] = "Y"
                    mappedNeighborL.append(activityD)
                else:
                    unmappedL.append(activityD)
                    activityD["neighbor_in_pdb"] = "N"
        #
        numLeft = maxActivity - len(mappedNeighborL)
        if numLeft > 0:
            unmappedL = sorted(unmappedL, key=lambda k: k["measurement_value"], reverse=True)
            retL = mappedNeighborL
            retL.extend(unmappedL[:numLeft])
            retL = sorted(retL, key=lambda k: k["measurement_value"], reverse=True)
        else:
            logger.debug("Mapped neighbor cofactors (%d) excluded unmapped (%d)", len(mappedNeighborL), len(unmappedL))
            retL = sorted(mappedNeighborL, key=lambda k: k["measurement_value"], reverse=True)

        return retL

    def __decodeComment(self, comment, separator="|"):
        dD = {}
        try:
            ti = iter(comment.split(separator))
            dD = {tup[1]: tup[0] for tup in zip(ti, ti)}
        except Exception:
            pass
        return dD
