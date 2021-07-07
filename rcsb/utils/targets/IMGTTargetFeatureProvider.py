##
#  File:           IMGTTargetFeatureProvider.py
#  Date:           7-Jul-2021 jdw
#
#  Updated:
#
##
"""
Accessors for IMGT (The International Immunogenetic Information System) target features
"""

import datetime
import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.targets.IMGTTargetProvider import IMGTTargetProvider

logger = logging.getLogger(__name__)


class IMGTTargetFeatureProvider(StashableBase):
    """Accessors for IMGT (The International Immunogenetic Information System) target features"""

    # Link out using the IMGT -
    # http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode=5w5m&Part=Chain
    #
    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirName = "IMGT-features"
        super(IMGTTargetFeatureProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__fD = self.__reload(self.__dirPath, useCache)
        #

    def testCache(self, minCount=20000):
        logger.info("IMGT feature count %d", len(self.__fD["features"]) if "features" in self.__fD else 0)
        if self.__fD and "features" in self.__fD and len(self.__fD["features"]) > minCount:
            return True
        else:
            return False

    def hasFeatures(self, rcsbInstanceId):
        """Return if features exist for the input instance identifier (auth_asym_id)

        Args:
            rcsbInstanceId (str): <pdbId (lower case)>.<auth_asym_id (case sensitive)>

        Returns:
            bool: True for success or False otherwise
        """
        return rcsbInstanceId in self.__fD["features"]

    def getFeatures(self, rcsbInstanceId):
        """Return features for the instance identifier (auth_asym_id)

        Args:
            rcsbInstanceId (str): <pdbId (lower case)>.<auth_asym_id (case sensitive)>

        Returns:
            list: list of feature dictionaries
        """
        try:
            return self.__fD["features"][rcsbInstanceId]
        except Exception:
            return None

    def __getFeatureDataPath(self):
        return os.path.join(self.__dirPath, "IMGT-feature-data.json")

    def reload(self):
        self.__fD = self.__reload(self.__dirPath, True)
        return True

    def __reload(self, dirPath, useCache):
        startTime = time.time()
        fD = {}
        featurePath = self.__getFeatureDataPath()
        #
        logger.info("useCache %r featurePath %r", useCache, featurePath)
        if useCache and self.__mU.exists(featurePath):
            fD = self.__mU.doImport(featurePath, fmt="json")
        else:
            fU = FileUtil()
            fU.mkdir(dirPath)
        # ---
        logger.info("Completed reload (useCache %r) at %s (%.4f seconds)", useCache, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return fD

    def buildFeatureList(self, useCache=True):
        """Build polymer instance feature list for IMGT annotations.

        Returns:
            bool: True for success or False otherwise

                    5w5m_B": {
                    "description": "FUSION-TNFRSF1B-GAMMA-1",
                    "domains": {
                        "C-DOMAIN|CH2|1": {
                        "geneAlleles": [
                            {
                                "taxName": "Homo sapiens",
                                "geneAllele": "IGHG4*01"
                            },
                            {
                                "taxName": "Homo sapiens",
                                "geneAllele": "IGHG4*03"
                            },
                            {
                                "taxName": "Homo sapiens",
                                "geneAllele": "IGHG4*04"
                            }
                        ],
                        "alignment": {
                            "begEntitySeqId": 7,
                            "endEntitySeqId": 116,
                            "begIMGTSeqId": "1",
                            "endIMGTSeqId": "105"
                        }
                        },
                        "C-DOMAIN|CH3|2": {
                        "geneAlleles": [
                            {
                                "taxName": "Homo sapiens",
                                "geneAllele": "IGHG4*01"
                            },
                            {
                                "taxName": "Homo sapiens",
                                "geneAllele": "IGHG4*04"
                            }
                        ],
                        "alignment": {
                            "begEntitySeqId": 117,
                            "endEntitySeqId": 221,
                            "begIMGTSeqId": "106",
                            "endIMGTSeqId": "209"
                        }
                        }
                    },
                    "proteinName": "IgG4 Sigma1 Fc",
                    "receptorType": "IG",
                    "receptorDescription": "FUSION-[TNFRSF1B]2-FC-GAMMA-1",
                    "species": "Homo sapiens (human)"
                },
        """
        rDL = []
        imgtP = IMGTTargetProvider(cachePath=self.__cachePath, useCache=useCache)
        #
        provenanceSource = "IMGT"
        refScheme = "PDB entity"
        assignVersion = imgtP.getVersion()
        #
        chainD = imgtP.getChains()
        #
        fTupL = [
            ("description", "IMGT_CHAIN_DESCRIPTION"),
            ("proteinName", "IMGT_PROTEIN_NAME"),
            ("receptorType", "IMGT_RECEPTOR_TYPE"),
            ("receptorDescription", "IMGT_RECEPTOR_DESCRIPTION"),
        ]
        ii = 1
        #
        for chainId, chD in chainD.items():
            entryId = chainId[:4]
            authAsymId = chainId.split("_")[1]
            # descriptive features -
            for fTup in fTupL:
                rD = {
                    "entry_id": entryId,
                    "auth_asym_id": authAsymId,
                    "type": fTup[1],
                    "feature_id": "IMGT_" + str(ii),
                    "name": chD[fTup[0]] if fTup[0] in chD else None,
                    "provenance_source": provenanceSource,
                    "reference_scheme": refScheme,
                    "assignment_version": assignVersion,
                    "feature_positions_beg_seq_id": None,
                    "feature_positions_end_seq_id": None,
                }
                rDL.append(rD)
                ii += 1
            # domain features -
            if "domains" not in chD:
                continue
            for domainId, dD in chD["domains"].items():
                dIdL = domainId.split("|")
                domainName = dIdL[0] + " " + dIdL[1]
                begSeqId = endSeqId = None
                if "alignment" in dD:
                    begSeqId = dD["alignment"]["begEntitySeqId"]
                    endSeqId = dD["alignment"]["endEntitySeqId"]
                else:
                    logger.debug("%r missing alignment for in %r", chainId, dD)
                #
                gaL = []
                if "geneAlleles" in dD:
                    for gD in dD["geneAlleles"]:
                        gaL.append(gD["geneAllele"])
                else:
                    logger.debug("%r missing gene and alleles for in %r", chainId, dD)
                #
                #
                rD = {
                    "entry_id": entryId,
                    "auth_asym_id": authAsymId,
                    "type": "IMGT_DOMAIN_NAME",
                    "feature_id": "IMGT_" + str(ii),
                    "name": domainName,
                    "provenance_source": provenanceSource,
                    "reference_scheme": refScheme,
                    "assignment_version": assignVersion,
                    "feature_positions_beg_seq_id": begSeqId,
                    "feature_positions_end_seq_id": endSeqId,
                }
                rDL.append(rD)
                ii += 1
                #
                for ga in gaL:
                    rD = {
                        "entry_id": entryId,
                        "auth_asym_id": authAsymId,
                        "type": "IMGT_GENE_ALLELE_NAME",
                        "feature_id": "IMGT_" + str(ii),
                        "name": ga,
                        "provenance_source": provenanceSource,
                        "reference_scheme": refScheme,
                        "assignment_version": assignVersion,
                        "feature_positions_beg_seq_id": begSeqId,
                        "feature_positions_end_seq_id": endSeqId,
                    }
                    rDL.append(rD)
                    ii += 1
        #
        qD = {}
        for rD in rDL:
            eId = rD["entry_id"] + "." + rD["auth_asym_id"]
            qD.setdefault(eId, []).append(rD)
        #
        logger.info("IMGT antibody chain features (%d)", len(qD))
        #
        fp = self.__getFeatureDataPath()
        tS = datetime.datetime.now().isoformat()
        vS = assignVersion
        ok = self.__mU.doExport(fp, {"version": vS, "created": tS, "features": qD}, fmt="json", indent=3)
        return ok
