##
#  File:           CARDTargetFeatureProvider.py
#  Date:           12-Jun-2021 jdw
#
#  Updated:
#
##
"""
Accessors and generators for CARD target feature data.
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


class CARDTargetFeatureProvider(StashableBase):
    """Accessors and generators for CARD target feature data."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirName = "CARD-features"
        super(CARDTargetFeatureProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fD = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=590):
        logger.info("CARD feature count %d", len(self.__fD["features"]) if "features" in self.__fD else 0)
        if self.__fD and "features" in self.__fD and len(self.__fD["features"]) > minCount:
            return True
        else:
            return False

    def hasFeatures(self, rcsbEntityId):
        return rcsbEntityId.upper() in self.__fD["features"]

    def getFeatures(self, rcsbEntityId):
        try:
            return self.__fD["features"][rcsbEntityId.upper()]
        except Exception:
            return None

    def reload(self):
        self.__fD = self.__reload(self.__dirPath, useCache=True)
        return True

    def __getFeatureDataPath(self):
        return os.path.join(self.__dirPath, "CARD-feature-data.json")

    def __reload(self, dirPath, **kwargs):
        startTime = time.time()
        fD = {}
        useCache = kwargs.get("useCache", True)
        ok = False
        featurePath = self.__getFeatureDataPath()
        #
        logger.info("useCache %r featurePath %r", useCache, featurePath)
        if useCache and self.__mU.exists(featurePath):
            fD = self.__mU.doImport(featurePath, fmt="json")
        else:
            fU = FileUtil()
            fU.mkdir(dirPath)
        # ---
        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return fD

    def buildFeatureList(self, sequenceMatchFilePath, useTaxonomy=False):
        """Build polymer entity feature list for the matching entities in the input sequence match file.

        Args:
            sequenceMatchFilePath (str): sequence match output file path
            useTaxonomy (bool): apply taxonomy criteria to filter matches. Defaults to False

        Returns:
            bool: True for success or False otherwise


        """
        rDL = []
        cardP = CARDTargetProvider(cachePath=self.__cachePath, useCache=False)
        mD = self.__mU.doImport(sequenceMatchFilePath, fmt="json")
        #
        provenanceSource = "CARD"
        refScheme = "PDB entity"
        assignVersion = cardP.getAssignmentVersion()
        for queryId, matchDL in mD.items():
            qCmtD = self.__decodeComment(queryId)
            modelId = qCmtD["modelId"]
            if not cardP.hasFeature(modelId):
                logger.info("Skipping CARD model %r", modelId)
                continue

            for matchD in matchDL:
                #
                if "alignedRegions" in matchD:
                    begSeqId = ",".join([str(arD["targetBegin"]) for arD in matchD["alignedRegions"]])
                    endSeqId = ",".join([str(arD["targetEnd"]) for arD in matchD["alignedRegions"]])
                else:
                    begSeqId = matchD["targetStart"]
                    endSeqId = matchD["targetEnd"]
                tCmtD = self.__decodeComment(matchD["target"])
                entryId = tCmtD["entityId"].split("_")[0]
                entityId = tCmtD["entityId"].split("_")[1]
                nm = cardP.getFeature(modelId, "modelName")
                descr = cardP.getFeature(modelId, "descr")
                featureId = cardP.getFeature(modelId, "id")
                rD = {
                    "entry_id": entryId,
                    "entity_id": entityId,
                    "type": "CARD_MODEL",
                    "feature_id": featureId,
                    "name": nm,
                    "description": descr,
                    "provenance_source": provenanceSource,
                    "reference_scheme": refScheme,
                    "assignment_version": assignVersion,
                    "feature_positions_beg_seq_id": begSeqId,
                    "feature_positions_end_seq_id": endSeqId,
                    "query_tax_name": matchD["queryTaxName"] if "queryTaxName" in matchD else None,
                    "target_tax_name": matchD["targetTaxName"] if "targetTaxName" in matchD else None,
                    "match_status": matchD["taxonomyMatchStatus"] if "taxonomyMatchStatus" in matchD else None,
                    "lca_tax_name": matchD["lcaTaxName"] if "lcaTaxName" in matchD else None,
                }
                rDL.append(rD)
        #
        qD = {}
        for rD in rDL:
            eId = rD["entry_id"] + "_" + rD["entity_id"]
            qD.setdefault(eId, []).append(rD)
        # --
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

        fp = self.__getFeatureDataPath()
        tS = datetime.datetime.now().isoformat()
        vS = datetime.datetime.now().strftime("%Y-%m-%d")
        ok = self.__mU.doExport(fp, {"version": vS, "created": tS, "features": fqD}, fmt="json", indent=3)
        return ok

    def __decodeComment(self, comment, separator="|"):
        dD = {}
        try:
            ti = iter(comment.split(separator))
            dD = {tup[1]: tup[0] for tup in zip(ti, ti)}
        except Exception:
            pass
        return dD
