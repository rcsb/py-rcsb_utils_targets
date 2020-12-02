##
#  File:           SAbDabTargetProvider.py
#  Date:           30-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for Thera-SAbDab(Therapeutic Structural Antibody Database) target data mappings
"""

import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class SAbDabTargetProvider(object):
    """Accessors for SAbDab target assignments."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "SAbDab-targets")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__oD = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=590):
        logger.info("SAbDab count %d", len(self.__oD))
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def __reload(self, dirPath, **kwargs):
        startTime = time.time()
        useCache = kwargs.get("useCache", True)
        targetUrl = kwargs.get("targetUrl", "http://opig.stats.ox.ac.uk/webapps/newsabdab/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv")
        ok = False
        fU = FileUtil()
        _, dumpFileName = os.path.split(targetUrl)
        #
        fU.mkdir(dirPath)
        dumpPath = os.path.join(dirPath, dumpFileName)
        fastaPath = os.path.join(dirPath, "SAbDab-targets.fasta")
        dataPath = os.path.join(dirPath, "SAbDab-data.json")
        #
        logger.info("useCache %r sabdabDumpPath %r", useCache, dumpPath)
        if useCache and self.__mU.exists(dataPath) and self.__mU.exists(fastaPath):
            oD = self.__mU.doImport(dataPath, fmt="json")
        else:
            logger.info("Fetching url %s path %s", targetUrl, dumpPath)
            ok = fU.get(targetUrl, dumpPath)
        #
        if useCache and fU.exists(fastaPath):
            pass
        else:
            oD = self.__convertDumpToFasta(dumpPath=dumpPath, fastaPath=fastaPath, dataPath=dataPath)
        # ---
        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return oD

    def __convertDumpToFasta(self, dumpPath="TheraSAbDab_SeqStruc_OnlineDownload.csv", fastaPath="antibody-seq.fasta", dataPath="antibody-data.json"):
        ok1 = ok2 = False
        oD = {}
        try:
            mU = MarshalUtil()
            rDL = mU.doImport(dumpPath, fmt="csv", rowFormat="dict")
            logger.debug("rD keys %r", list(rDL[0].keys()))
            seqObj = {}
            for rD in rDL:
                tS = "|".join(rD[ky].strip() for ky in ["Therapeutic", "Format", "CH1 Isotype", "VD LC", "Highest_Clin_Trial (Jan '20)", "Est. Status"])
                tS = "".join(tS.split())
                #
                oD[rD["Therapeutic"]] = {
                    kTup[1]: rD[kTup[0]]
                    for kTup in [
                        ("Therapeutic", "name"),
                        ("Format", "format"),
                        ("CH1 Isotype", "chiIsotype"),
                        ("VD LC", "VD_LC"),
                        ("Highest_Clin_Trial (Jan '20)", "maxPhase"),
                        ("Est. Status", "status"),
                    ]
                }
                hSeq = rD["Heavy Sequence"] if rD["Heavy Sequence"] != "na" else None
                lSeq = rD["Light Sequence"] if rD["Light Sequence"] != "na" else None
                if hSeq:
                    tS = rD["Therapeutic"] + "|heavy"
                    seqObj[tS] = {"sequence": hSeq.strip(), "therapeutic": rD["Therapeutic"], "chain": "heavy"}
                if lSeq:
                    tS = rD["Therapeutic"] + "|light"
                    seqObj[tS] = {"sequence": lSeq.strip(), "therapeutic": rD["Therapeutic"], "chain": "light"}
                #
            ok1 = mU.doExport(fastaPath, seqObj, fmt="fasta", makeComment=True)
            ok2 = mU.doExport(dataPath, oD, fmt="json", indent=3)
            logger.info("Exporting SAbDab %d sequences in %r status %r", len(seqObj), fastaPath, ok1 and ok2)
        except Exception as e:
            logger.exception("Failing for %r with %s", dumpPath, str(e))
        return oD