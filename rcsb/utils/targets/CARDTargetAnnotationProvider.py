##
#  File:           CARDTargetAnnotationProvider.py
#  Date:           6-Mar-2023 dwp
#
#  Updated:
#
##
"""
Accessors and generators for CARD target annotation data.
"""

import datetime
import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.targets.CARDTargetProvider import CARDTargetProvider

logger = logging.getLogger(__name__)


class CARDTargetAnnotationProvider(StashableBase):
    """Accessors and generators for CARD target annotation data."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirName = "CARD-annotations"
        super(CARDTargetAnnotationProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__cardP = CARDTargetProvider(cachePath=self.__cachePath, useCache=False)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fD = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=590):
        logger.info("CARD annotation count %d", len(self.__fD["annotations"]) if "annotations" in self.__fD else 0)
        if self.__fD and "annotations" in self.__fD and len(self.__fD["annotations"]) > minCount:
            return True
        else:
            return False

    def hasAnnotation(self, rcsbEntityId):
        return rcsbEntityId.upper() in self.__fD["annotations"]

    def getAnnotation(self, rcsbEntityId):
        try:
            return self.__fD["annotations"][rcsbEntityId.upper()]
        except Exception:
            pass
        return None

    def reload(self):
        self.__fD = self.__reload(self.__dirPath, useCache=True)
        return True

    def __getAnnotationDataPath(self):
        return os.path.join(self.__dirPath, "CARD-annotation-data.json")

    def __reload(self, dirPath, **kwargs):
        startTime = time.time()
        fD = {}
        useCache = kwargs.get("useCache", True)
        ok = False
        annotationPath = self.__getAnnotationDataPath()
        #
        logger.info("useCache %r annotationPath %r", useCache, annotationPath)
        if useCache and self.__mU.exists(annotationPath):
            fD = self.__mU.doImport(annotationPath, fmt="json")
        else:
            fU = FileUtil()
            fU.mkdir(dirPath)
        # ---
        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return fD

    def buildAnnotationList(self, sequenceMatchFilePath, useTaxonomy=False):
        """Build polymer entity annotation list for the matching entities in the input sequence match file.

        Args:
            sequenceMatchFilePath (str): sequence match output file path
            useTaxonomy (bool): apply taxonomy criteria to filter matches. Defaults to False

        Returns:
            bool: True for success or False otherwise
        """
        rDL = []
        mD = self.__mU.doImport(sequenceMatchFilePath, fmt="json")
        #
        provenanceSource = "matching CARD Protein Homolog Models (PHM)"
        refScheme = "PDB entity"
        assignVersion = self.__cardP.getAssignmentVersion()
        for queryId, matchDL in mD.items():
            qCmtD = self.__decodeComment(queryId)
            modelId = qCmtD["modelId"]
            if not self.__cardP.hasModel(modelId):
                logger.info("Skipping CARD model %r", modelId)
                continue

            for matchD in matchDL:
                #
                fpL = []
                tCmtD = self.__decodeComment(matchD["target"])
                entryId = tCmtD["entityId"].split("_")[0]
                entityId = tCmtD["entityId"].split("_")[1]
                seqIdPct = matchD["sequenceIdentity"]
                bitScore = matchD["bitScore"]
                nm = self.__cardP.getModelKey(modelId, "modelName")
                descr = self.__cardP.getModelKey(modelId, "descr")
                cvTermId = self.__cardP.getModelKey(modelId, "cvTermId")
                annotationId = self.__cardP.getModelKey(modelId, "accession")
                familyCvTermId = self.__cardP.getModelKey(modelId, "familyCvTermId")
                familyName = self.__cardP.getModelKey(modelId, "familyName")
                familyAnnotationId = self.__cardP.getModelKey(modelId, "familyAccession")
                familyDescription = self.__cardP.getModelKey(modelId, "familyDescription")
                drugClasses = self.__cardP.getModelKey(modelId, "drugClasses")
                resistanceMechanism = self.__cardP.getModelKey(modelId, "resistanceMechanism")
                rD = {
                    "entry_id": entryId,
                    "entity_id": entityId,
                    "seq_id_pct": seqIdPct,
                    "bit_score": bitScore,
                    # "type": "CARD",  # Assigned in DictMethodEntityHelper
                    "annotation_id": "ARO:" + annotationId,
                    "card_aro_cvterm_id": cvTermId,
                    "name": nm,
                    "family_annotation_id": "ARO:" + familyAnnotationId,
                    "family_card_aro_cvterm_id": familyCvTermId,
                    "family_name": familyName,
                    "family_description": familyDescription,
                    "card_aro_drug_class": drugClasses,
                    "card_aro_resistance_mechanism": resistanceMechanism,
                    "description": descr,
                    "provenance_source": provenanceSource,
                    "reference_scheme": refScheme,
                    "assignment_version": assignVersion,
                    "annotation_positions": fpL,
                    "query_tax_name": matchD["queryTaxName"] if "queryTaxName" in matchD else None,
                    "target_tax_name": matchD["targetTaxName"] if "targetTaxName" in matchD else None,
                    "match_status": matchD["taxonomyMatchStatus"] if "taxonomyMatchStatus" in matchD else None,
                    "lca_tax_name": matchD["lcaTaxName"] if "lcaTaxName" in matchD else None,
                }
                rDL.append(rD)
        #
        qD = {}
        dD = {}
        for rD in rDL:
            eId = rD["entry_id"] + "_" + rD["entity_id"]
            fId = rD["annotation_id"]
            if (eId, fId) in dD:
                continue
            dD[(eId, fId)] = True
            qD.setdefault(eId, []).append(rD)
        # --
        # If useTaxonomy filter is True, resulting dictionary contains ~100 fewer perfect matches (as of March 2023)
        if useTaxonomy:
            fqD = {}
            for eId, rDL in qD.items():
                mL = []
                oL = []
                for rD in rDL:
                    if "match_status" not in rD:
                        continue
                    if rD["match_status"] in ["matched", "alternate strain", "alternate strain/serotype/isolate/genotype", "query is ancestor", "query is descendant"]:
                        mL.append(rD)
                    elif rD["match_status"] == "orthologous match (by lca)":
                        oL.append(rD)
                    else:
                        continue
                if not mL and not oL:
                    continue
                logger.debug("eId %r mL (%d) oL (%d)", eId, len(mL), len(oL))
                fqD.setdefault(eId, []).extend(mL if mL else oL)
        else:
            fqD = qD
        # --
        # Determine whether entities demonstrate a perfect match or match multiple CARD annotations
        cqD = {}
        for eId, rDL in fqD.items():
            if len(rDL) == 1:
                rD = rDL[0]
                rD["perfect_match"] = "Y"
                cqD[eId] = rD
            else:
                # Sort list of matches by highest sequence identity percent and bit score
                rDLSorted = sorted(rDL, key=lambda k: (-k["seq_id_pct"], -k["bit_score"]))
                rD = rDLSorted[0]
                # Check for rD that has 100% match, and if so, set this a perfect_match
                if rD["seq_id_pct"] == 100.0:
                    rD["perfect_match"] = "Y"
                else:
                    rD["perfect_match"] = "N"
                rDL[0]["perfect_match"] = "N"  # indicates there were multiple CARDs matching the given entity
                cqD[eId] = rD
        # --
        fp = self.__getAnnotationDataPath()
        tS = datetime.datetime.now().isoformat()
        vS = datetime.datetime.now().strftime("%Y-%m-%d")
        ok = self.__mU.doExport(fp, {"version": vS, "created": tS, "taxonomyFilter": useTaxonomy, "annotations": cqD}, fmt="json", indent=3)
        return ok

    def __decodeComment(self, comment, separator="|"):
        dD = {}
        try:
            ti = iter(comment.split(separator))
            dD = {tup[1]: tup[0] for tup in zip(ti, ti)}
        except Exception:
            pass
        return dD

    def getLineage(self, aroId):
        """Return the lineage (all parents + the requested aroId) of the given aroId.

        Args:
            aroId (str): ARO ID in the form of "ARO:3001059"

        Returns:
            list: list of dictionaries containing the "id" and "name" of each parent
        """
        return self.__cardP.getLineage(aroId)
