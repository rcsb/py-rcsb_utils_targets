##
#  File:           CARDTargetProvider.py
#  Date:           27-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for CARD target assignments.

"""

import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class CARDTargetProvider:
    """Accessors for CARD target assignments."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "CARD-targets")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__oD, self.__version = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=3000):
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def __reload(self, dirPath, **kwargs):
        oD = None
        version = None
        startTime = time.time()
        useCache = kwargs.get("useCache", True)
        #
        # CARDDumpUrl = kwargs.get("CARDDumpUrl", "https://card.mcmaster.ca/latest/data/broadstreet-v3.1.0.tar.bz2")
        cardDumpUrl = kwargs.get("CARDDumpUrl", "https://card.mcmaster.ca/latest/data")
        ok = False
        fU = FileUtil()
        cardDumpFileName = "card-data.tar.bz2"
        cardDumpPath = os.path.join(dirPath, cardDumpFileName)
        cardDumpDirPath = os.path.join(dirPath, "dump")
        #
        fU.mkdir(dirPath)
        cardDataPath = os.path.join(dirPath, "card-select-data.json")
        cardFastaPath = os.path.join(dirPath, "card-targets.fasta")
        #
        logger.info("useCache %r CARDDumpPath %r", useCache, cardDumpPath)
        if useCache and self.__mU.exists(cardDataPath):
            qD = self.__mU.doImport(cardDataPath, fmt="json")
            version = qD["version"]
            oD = qD["data"]
        else:
            logger.info("Fetching url %s path %s", cardDumpUrl, cardDumpPath)
            ok = fU.get(cardDumpUrl, cardDumpPath)
            fU.mkdir(cardDumpDirPath)
            fU.uncompress(cardDumpPath, outputDir=cardDumpDirPath)
            fU.unbundleTarfile(os.path.join(cardDumpDirPath, cardDumpFileName[:-4]), dirPath=cardDumpDirPath)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            oD, version = self.__parseCardData(os.path.join(cardDumpDirPath, "card.json"))
            ok1 = self.__mU.doExport(cardDataPath, {"version": version, "data": oD}, fmt="json", indent=3)
            ok2 = self.__exportCardFasta(cardFastaPath, oD)
            ok = ok1 and ok2
            logger.info("Export CARD data (%d) status %r", len(oD), ok)
        # ---
        return oD, version

    def __exportCardFasta(self, fastaPath, cardD):
        """Export a CARD sequence target fasta file

        Args:
            fastaPath (str): fasta output file path
            cardD (dict): card selected data dictionary

        Returns:
            (bool): True for success or False otherwise
        """
        sD = {}
        try:
            for modelId, tD in cardD.items():
                # aroAcc = tD["accession"]
                aroId = tD["id"]
                if "sequences" not in tD:
                    continue
                for qD in tD["sequences"]:
                    sId = qD["seqId"]
                    seq = qD["sequence"]
                    taxId = qD["taxId"]
                    sD[sId] = {"sequence": seq, "modelId": modelId, "aroId": aroId, "seqId": sId, "taxId": taxId}

            ok = self.__mU.doExport(fastaPath, sD, fmt="fasta", makeComment=True)
            logger.info("Export CARD fasta (%d) status %r", len(sD), ok)
        except Exception as e:
            logger.exception("Failing for model %r tD %r with %s", modelId, tD, str(e))
        return ok

    def __parseCardData(self, filePath):
        """Parse CARD target data

        Args:
            filePath (str): card json data file

        Returns:
            (dict, string): card selected data dictionary, card version string

        """
        try:
            oD = {}
            version = None
            cD = self.__mU.doImport(filePath, fmt="json")
            logger.info("CARD model count (%d)", len(cD))
            for modelId, mD in cD.items():
                if modelId.startswith("_"):
                    version = mD
                    continue
                oD[modelId] = {}
                for kTup in [
                    ("ARO_accession", "accession"),
                    ("ARO_id", "id"),
                    ("ARO_name", "name"),
                    ("ARO_description", "descr"),
                    ("model_name", "modelName"),
                    ("model_type", "modelType"),
                ]:
                    if kTup[0] in mD:
                        oD[modelId][kTup[1]] = mD[kTup[0]]

                try:
                    if "model_sequences" in mD:
                        for seqId, tD in mD["model_sequences"]["sequence"].items():
                            oD[modelId].setdefault("sequences", []).append(
                                {"seqId": seqId, "sequence": tD["protein_sequence"]["sequence"], "taxId": tD["NCBI_taxonomy"]["NCBI_taxonomy_id"]}
                            )
                except Exception as e:
                    logger.exception("Failing with %s", str(e))
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD, version