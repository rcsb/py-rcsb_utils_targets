##
#  File:           SAbDabTargetProvider.py
#  Date:           18-Jun-2021 jdw
#
#  Updated:
#
##
"""
Accessors for Thera-SAbDab(Therapeutic Structural Antibody Database) target data.
"""

import datetime
import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class SAbDabTargetProvider(object):
    """Accessors for SAbDab and Thera-SAbDab(Therapeutic Structural Antibody Database) target data.

       See: Dunbar, J., Krawczyk, K. et al (2014). Nucleic Acids Res. 42. D1140-D1146
    """

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "SAbDab-features")
        #
        self.__assignVersion = None
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__oD, self.__aD, self.__dumpPath, self.__assignVersion = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=590):
        logger.info("SAbDab count %d", len(self.__oD))
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def getFeatures(self, therapeuticName, featureKey):
        fL = []
        try:
            fS = self.__oD[therapeuticName][featureKey]
            if ";" in fS:
                fL = fS.split(";")
            else:
                fL = [fS]
        except Exception:
            fL = []
        return fL

    def getAssignment(self, instanceId, featureKey):
        """Return the value of the key feature for the input instance identifier.

        Args:
            instanceId (str): instance identifier '<pdbId>|<authAsymId>'
            featureKey (str): assignment feature key: pdb|Hchain|Lchain|model|antigen_chain|antigen_type|
                              antigen_het_name|antigen_name|heavy_subclass|light_subclass|light_ctype)

        Returns:
            str:  feature value or None
        """
        fVal = None
        try:
            fVal = self.__aD[instanceId][featureKey]
        except Exception:
            fVal = None
        return fVal

    def getAssignmentVersion(self):
        return self.__assignVersion

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
        dataPath = os.path.join(dirPath, "sabdab-data.json")
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
            tD = {}
            for rD in rDL:
                tD[rD["Therapeutic"]] = {
                    kTup[1]: rD[kTup[0]] if rD[kTup[0]] not in ["na", "na;na"] else None
                    for kTup in [
                        ("Therapeutic", "antibodyName"),
                        ("Format", "antiBodyFormat"),
                        ("CH1 Isotype", "ch1Isotype"),
                        ("VD LC", "VD_LC"),
                        ("Highest_Clin_Trial (Jan '20)", "maxClinicalPhase"),
                        ("Est. Status", "status"),
                        ("Target", "target"),
                        ("Conditions Approved", "conditionsApproved"),
                        ("Conditions Active", "conditionsActive"),
                    ]
                }
            aD = self.__reloadAssignments(dirPath, **kwargs)
            #
            tS = datetime.datetime.now().isoformat()
            vS = datetime.datetime.now().strftime("%Y-%m-%d")
            oD = {"version": vS, "created": tS, "identifiers": tD, "assignments": aD}
            ok = self.__mU.doExport(dataPath, oD, fmt="json", indent=3)
            logger.info("Exporting (%d) Thera-SAbDab data records and (%d) SAbDab assignments in %r status %r", len(oD["identifiers"]), len(oD["assignments"]), dataPath, ok)

        # ---
        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return oD["identifiers"], oD["assignments"], dumpPath, oD["version"]

    def __reloadAssignments(self, dirPath, **kwargs):
        """Fetch and read

        Args:
            dirPath ([type]): [description]

        Returns:
            [type]: [description]

        """
        startTime = time.time()
        aD = {}
        try:
            targetUrl = kwargs.get("assignmentUrl", "http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/summary/all")
            fU = FileUtil()
            dumpFileName = "sabdab_summary_all.tsv"
            #
            fU.mkdir(dirPath)
            dumpPath = os.path.join(dirPath, dumpFileName)
            logger.info("Fetching url %s path %s", targetUrl, dumpPath)
            ok = fU.get(targetUrl, dumpPath)
            rDL = self.__mU.doImport(dumpPath, fmt="tdd", rowFormat="dict")
            logger.info("SAbDab raw records (%d)", len(rDL))
            logger.debug("rD keys %r", list(rDL[0].keys()))
            kyL = ["pdb", "Hchain", "Lchain", "model", "antigen_chain", "antigen_type", "antigen_het_name", "antigen_name", "heavy_subclass", "light_subclass", "light_ctype"]

            for rD in rDL:
                pdbId = rD["pdb"] if rD["pdb"] and rD["pdb"] != "NA" else None
                authAsymIdH = rD["Hchain"] if rD["Hchain"] and rD["Hchain"] != "NA" else None
                authAsymIdL = rD["Lchain"] if rD["Lchain"] and rD["Lchain"] != "NA" else None
                if pdbId and authAsymIdH:
                    aD[pdbId + "|" + authAsymIdH] = {k: v for k, v in rD.items() if v and v != "NA" and k in kyL}
                if pdbId and authAsymIdL:
                    aD[pdbId + "|" + authAsymIdL] = {k: v for k, v in rD.items() if v and v != "NA" and k in kyL}

            logger.info("Fetched (%d) SAbDab assignment records.", len(aD))
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        logger.info("Completed reload (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return aD

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
                if hSeq:
                    cD = {"sequence": hSeq.strip(), "therapeutic": rD["Therapeutic"], "chain": "heavy"}
                seqId = ""
                cL = []
                for k, v in cD.items():
                    if k in ["sequence"]:
                        continue
                    cL.append(str(v))
                    cL.append(str(k))
                #
                seqId = "|".join(cL)
                sD[seqId] = cD
            for rD in rDL:
                lSeq = rD["Light Sequence"] if rD["Light Sequence"] != "na" else None
                if lSeq:
                    cD = {"sequence": lSeq.strip(), "therapeutic": rD["Therapeutic"], "chain": "light"}
                seqId = ""
                cL = []
                for k, v in cD.items():
                    if k in ["sequence"]:
                        continue
                    cL.append(str(v))
                    cL.append(str(k))
                #
                seqId = "|".join(cL)
                sD[seqId] = cD
                #
            ok = self.__mU.doExport(fastaPath, sD, fmt="fasta", makeComment=True)
            logger.info("Exporting SAbDab %d fasta sequences in %r status %r", len(sD), fastaPath, ok)
        except Exception as e:
            logger.exception("Failing for %r with %s", dumpPath, str(e))
        return ok
