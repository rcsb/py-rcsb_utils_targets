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

import datetime
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
        useCache = kwargs.get("useCache", True)
        self.__dirName = "SAbDab-features"
        super(SAbDabTargetFeatureProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fD = self.__reload(self.__dirPath, useCache)
        #

    def testCache(self, minCount=500):
        logger.info(
            "Therapeutic SAbDab feature count %d Assignment count %d",
            len(self.__fD["features"]) if "features" in self.__fD else 0,
            len(self.__fD["assignments"]) if "assignments" in self.__fD else 0,
        )
        if self.__fD and "features" in self.__fD and len(self.__fD["features"]) > minCount and self.__fD and "assignments" in self.__fD and len(self.__fD["assignments"]):
            return True
        else:
            return False

    def getVersion(self):
        try:
            return self.__fD["version"]
        except Exception:
            pass
        return None

    def hasFeatures(self, rcsbEntityId):
        return rcsbEntityId.upper() in self.__fD["features"]

    def getFeatures(self, rcsbEntityId):
        try:
            return self.__fD["features"][rcsbEntityId.upper()]
        except Exception:
            return []

    def getAssignment(self, instanceId, featureKey):
        """Return the value of the key feature for the input instance identifier.

        Args:
            instanceId (str): instance identifier '<pdbId>.<authAsymId>'
            featureKey (str): assignment feature key: one of pdb|Hchain|Lchain|model|antigen_chain|antigen_type|
                              antigen_het_name|antigen_name|heavy_subclass|light_subclass|light_ctype)

        Returns:
            str:  feature value or None
        """
        fVal = None
        try:
            fVal = self.__fD["assignments"][instanceId][featureKey]
        except Exception:
            fVal = None
        return fVal

    def hasAssignment(self, instanceId):
        """Return if assignment data is available for the input instance.

        Args:
            instanceId (str): instance identifier '<pdbId>.<authAsymId>'

        Returns:
            bool: True for success or False otherwise
        """
        return instanceId in self.__fD["assignments"]

    def __getFeatureDataPath(self):
        return os.path.join(self.__dirPath, "sabdab-feature-data.json")

    def reload(self):
        self.__fD = self.__reload(self.__dirPath, True)
        return True

    def __reload(self, dirPath, useCache):
        startTime = time.time()
        fD = {}

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
        #
        # - sort out if we match light and heavy chains
        #
        iD = {}
        fullMatchD = {}
        for queryId, matchDL in mD.items():
            qCmtD = self.__decodeComment(queryId)
            # Tanezumab|therapeutic|light|chain
            thName = qCmtD["therapeutic"]
            chainType = qCmtD["chain"]
            for matchD in matchDL:
                tCmtD = self.__decodeComment(matchD["target"])
                entryId = tCmtD["entityId"].split("_")[0]
                entityId = tCmtD["entityId"].split("_")[1]
                iD[(thName, chainType, entryId)] = entityId
        logger.info("Match index length (%d)", len(iD))
        for (thName, chainType, entryId), entityId in iD.items():
            if chainType == "light":
                continue
            if (thName, "light", entryId) in iD:
                fullMatchD[(thName, "heavy", entryId, entityId)] = True
                lEntityId = iD[(thName, "light", entryId)]
                fullMatchD[(thName, "light", entryId, lEntityId)] = True
        logger.info("Antibody entity match length (%d)", len(fullMatchD))
        #
        # - Add features for full matches -
        for queryId, matchDL in mD.items():
            qCmtD = self.__decodeComment(queryId)
            # Tanezumab|therapeutic|light|chain
            thName = qCmtD["therapeutic"]
            chainType = qCmtD["chain"]
            #
            for matchD in matchDL:
                fpL = []
                if "alignedRegions" in matchD:
                    fpL = [{"beg_seq_id": arD["targetBegin"], "end_seq_id": arD["targetEnd"]} for arD in matchD["alignedRegions"]]
                else:
                    fpL = [{"beg_seq_id": matchD["targetBegin"], "end_seq_id": matchD["targetEnd"]}]
                #
                tCmtD = self.__decodeComment(matchD["target"])
                entryId = tCmtD["entityId"].split("_")[0]
                entityId = tCmtD["entityId"].split("_")[1]
                if (thName, chainType, entryId, entityId) not in fullMatchD:
                    continue
                ii = 1
                for fType, fKy in [
                    ("SABDAB_ANTIBODY_NAME", "antibodyName"),
                    ("SABDAB_ANTIBODY_FORMAT", "antibodyFormat"),
                    ("SABDAB_ANTIBODY_CH1_ISOTYPE", "ch1IsoType"),
                    ("SABDAB_ANTIBODY_LIGHT_CHAIN_TYPE", "VD_LC"),
                    ("SABDAB_ANTIBODY_TARGET", "target"),
                ]:
                    if fType == "Antibody_Light_Chain_Type" and chainType == "heavy":
                        continue
                    fVL = stP.getFeatures(thName, fKy)
                    if not fVL:
                        continue
                    for fV in fVL:
                        rD = {
                            "entry_id": entryId,
                            "entity_id": entityId,
                            "type": fType,
                            "feature_id": thName + "_" + chainType + "_" + str(ii),
                            "name": fV,
                            "provenance_source": provenanceSource,
                            "reference_scheme": refScheme,
                            "assignment_version": assignVersion,
                            "feature_positions": fpL,
                        }
                        rDL.append(rD)
                        ii += 1
        #
        qD = {}
        for rD in rDL:
            eId = rD["entry_id"] + "_" + rD["entity_id"]
            qD.setdefault(eId, []).append(rD)
        #
        logger.info("Antibody matches (%d)", len(qD))
        #
        fp = self.__getFeatureDataPath()
        tS = datetime.datetime.now().isoformat()
        vS = datetime.datetime.now().strftime("%Y-%m-%d")
        ok = self.__mU.doExport(fp, {"version": vS, "created": tS, "features": qD, "assignments": stP.getAssignments()}, fmt="json", indent=3)
        return ok

    def __decodeComment(self, comment, separator="|"):
        dD = {}
        try:
            ti = iter(comment.split(separator))
            dD = {tup[1]: tup[0] for tup in zip(ti, ti)}
        except Exception:
            pass
        return dD
