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
        self.__oD, self.__dumpPath = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=590):
        logger.info("SAbDab count %d", len(self.__oD))
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def __reload(self, dirPath, **kwargs):
        startTime = time.time()
        oD = {}
        useCache = kwargs.get("useCache", True)
        targetUrl = kwargs.get("targetUrl", "http://opig.stats.ox.ac.uk/webapps/newsabdab/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv")
        #
        ok = False
        fU = FileUtil()
        _, dumpFileName = os.path.split(targetUrl)
        #
        fU.mkdir(dirPath)
        dumpPath = os.path.join(dirPath, dumpFileName)
        dataPath = os.path.join(dirPath, "SAbDab-data.json")
        #
        logger.info("useCache %r sabdabDumpPath %r", useCache, dumpPath)
        if useCache and self.__mU.exists(dataPath):
            oD = self.__mU.doImport(dataPath, fmt="json")
        else:
            logger.info("Fetching url %s path %s", targetUrl, dumpPath)
            ok = fU.get(targetUrl, dumpPath)
            #
            rDL = self.__mU.doImport(dumpPath, fmt="csv", rowFormat="dict")
            logger.debug("rD keys %r", list(rDL[0].keys()))
            for rD in rDL:
                oD[rD["Therapeutic"]] = {
                    kTup[1]: rD[kTup[0]] if rD[kTup[0]] not in ["na", "na;na"] else None
                    for kTup in [
                        ("Therapeutic", "antibodyName"),
                        ("Format", "antiBodyFormat"),
                        ("CH1 Isotype", "ch1Isotype"),
                        ("VD LC", "VD_LC"),
                        ("Highest_Clin_Trial (Jan '20)", "maxClinicalPhase"),
                        ("Est. Status", "status"),
                    ]
                }

            ok = self.__mU.doExport(dataPath, oD, fmt="json", indent=3)
            logger.info("Exporting SAbDab %d data records in %r status %r", len(oD), dataPath, ok)

        # ---
        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return oD, dumpPath

    def exportFasta(self, fastaPath):
        ok = self.__convertDumpToFasta(dumpPath=self.__dumpPath, fastaPath=fastaPath)
        return ok

    #
    def __convertDumpToFasta(self, dumpPath, fastaPath):
        ok = False
        try:
            rDL = self.__mU.doImport(dumpPath, fmt="csv", rowFormat="dict")
            logger.debug("rD keys %r", list(rDL[0].keys()))
            sD = {}
            for rD in rDL:
                hSeq = rD["Heavy Sequence"] if rD["Heavy Sequence"] != "na" else None
                lSeq = rD["Light Sequence"] if rD["Light Sequence"] != "na" else None
                if hSeq:
                    cD = {"sequence": hSeq.strip(), "therapeutic": rD["Therapeutic"], "chain": "heavy"}
                if lSeq:
                    cD = {"sequence": lSeq.strip(), "therapeutic": rD["Therapeutic"], "chain": "light"}
                seqId = ""
                cL = []
                for k, v in cD.items():
                    if k in ["sequence"]:
                        continue
                    cL.append(str(v))
                    cL.append(str(k))
                seqId = "|".join(cL)
                sD[seqId] = cD
                #
            ok = self.__mU.doExport(fastaPath, sD, fmt="fasta", makeComment=True)
            logger.info("Exporting SAbDab %d fasta sequences in %r status %r", len(sD), fastaPath, ok)
        except Exception as e:
            logger.exception("Failing for %r with %s", dumpPath, str(e))
        return ok
