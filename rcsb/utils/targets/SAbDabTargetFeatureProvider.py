##
#  File:           SAbDabTargetFeatureProvider.py
#  Date:           11-Jun-2021 jdw
#
#  Updated:
#
##
"""
Accessors for Thera-SAbDab(Therapeutic Structural Antibody Database) target features.
"""

import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.targets.SAbDabTargetProvider import SAbDabTargetProvider

logger = logging.getLogger(__name__)


class SAbDabTargetFeatureProvider(StashableBase):
    """Accessors for Thera-SAbDab(Therapeutic Structural Antibody Database) target features."""

    # Link out using the INN therapeutic name -
    # http://opig.stats.ox.ac.uk/webapps/newsabdab/therasabdab/search/?therapeutic=Coltuximab
    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirName = "SAbDab-features"
        super(SAbDabTargetFeatureProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fD = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=590):
        logger.info("SAbDab feature count %d", len(self.__fD["features"]) if "features" in self.__fD else 0)
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

    def __getFeatureDataPath(self):
        return os.path.join(self.__dirPath, "SAbDab-feature-data.json")

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

    def buildFeatureList(self, sequenceMatchFilePath):
        """Build polymer entity feature list for the matching entities in the input sequence match file.

        Args:
            sequenceMatchFilePath (str): sequence match output file path

        Returns:
            bool: True for success or False otherwise

        """
        rDL = []
        stP = SAbDabTargetProvider(cachePath=self.__cachePath, useCache=False)
        mD = self.__mU.doImport(sequenceMatchFilePath, fmt="json")
        #
        provenanceSource = "SAbDab"
        refScheme = "PDB entity"
        assignVersion = stP.getAssignmentVersion()
        for queryId, matchDL in mD.items():
            thName = queryId.split("|")[0]
            for matchD in matchDL:
                begSeqId = matchD["targetStart"]
                endSeqId = matchD["targetEnd"]
                tL = matchD["target"].split("|")
                entryId = tL[0].split("_")[0]
                entityId = tL[0].split("_")[1]
                for fType, fKy in [
                    ("Antibody_Name", "antibodyName"),
                    ("Antibody_Format", "antibodyFormat"),
                    ("Antibody_CH1_IsoType", "ch1IsoType"),
                    ("Antibody_Light_Chain_Type", "VD_LC"),
                    ("Antibody_Target", "target"),
                ]:
                    fVL = stP.getFeatures(thName, fKy)
                    if not fVL:
                        continue
                    for ii, fV in enumerate(fVL):
                        rD = {
                            "entry_id": entryId,
                            "entity_id": entityId,
                            "type": fType,
                            "feature_id": "SAbDab" + "_" + str(ii + 1),
                            "name": fV,
                            "provenance_source": provenanceSource,
                            "reference_scheme": refScheme,
                            "assignment_version": assignVersion,
                            "feature_positions_beg_seq_id": begSeqId,
                            "feature_positions_end_seq_id": endSeqId,
                        }
                        rDL.append(rD)
        #
        qD = {}
        for rD in rDL:
            eId = rD["entry_id"] + "_" + rD["entity_id"]
            qD.setdefault(eId, []).append(rD)
        fp = self.__getFeatureDataPath()
        tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        vS = time.strftime("%Y-%m-%d", time.localtime())
        ok = self.__mU.doExport(fp, {"version": vS, "created": tS, "features": qD}, fmt="json", indent=3)
        return ok
