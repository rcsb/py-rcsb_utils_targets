##
#  File:           SAbDabTargetProvider.py
#  Date:           18-Jun-2021 jdw
#
#  Updated:
#   9-Feb-2023 aae  Find Highest_Clin_Trial column regardless of month
#   1-Jul-2024 dwp  Update SAbDab data parsing (following change in source data headers)
#  22-Jul-2024 dwp  Re-update SAbDab data parsing (following change in source data headers)
#   2-Jul-2026 dwp  Use backup SAbDab summary file and TheraSAbDab file during transition phase to SAbDab2
#  19-Oct-2026 dwp  Switch SAbDab assignments to the SAbDab2 JSON API (rcsb-pdb-annotations)
##
"""
Accessors for Thera-SAbDab(Therapeutic Structural Antibody Database) target data.
"""

import datetime
import logging
import os.path
import re
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
        logger.info("Therapeutic SAbDab count %d", len(self.__oD))
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
            instanceId (str): instance identifier '<pdbId>.<authAsymId>'
            featureKey (str): assignment feature key: one of pdb|model|antigen_name|heavy_subclass|light_subclass|light_ctype
                              (antigen_name is ' | ' joined; light_ctype is Kappa|Lambda)

        Returns:
            str:  feature value or None
        """
        fVal = None
        try:
            fVal = self.__aD[instanceId][featureKey]
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
        return instanceId in self.__aD

    def getAssignments(self):
        return self.__aD

    def getAssignmentVersion(self):
        return self.__assignVersion

    def reload(self):
        self.__oD, self.__aD, self.__dumpPath, self.__assignVersion = self.__reload(self.__dirPath, useCache=True)

    def __reload(self, dirPath, **kwargs):
        startTime = time.time()
        oD = {}
        useCache = kwargs.get("useCache", True)
        # The old newsabdab Thera-SAbDab URL now serves an HTML page; this is the last known location of the CSV
        targetUrl = kwargs.get("targetUrl", "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv")
        targetFallbackUrl = kwargs.get(
            "targetFallbackUrl",
            "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets_stash/refs/heads/development/stash/SAbDab-backup/TheraSAbDab_SeqStruc_OnlineDownload.csv",
        )
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
            if not ok:
                logger.error("Fetching failed for Thera-SAbDab target data. Trying fallback %s", targetFallbackUrl)
                ok = fU.get(targetFallbackUrl, dumpPath)
                if not ok:
                    raise ValueError("Fetching failed for fallback Thera-SAbDab target data")
            #
            rDL = self.__mU.doImport(dumpPath, fmt="csv", rowFormat="dict")
            logger.debug("rD keys %r", list(rDL[0].keys()))
            tD = {}
            for rD in rDL:
                qD = {}
                # Highest_Clin_Trial column name changes each month
                clinTrialCol = next((k for k in list(rD.keys()) if "Highest_Clin_Trial" in k), "Highest_Clin_Trial")
                #
                for kTup in [
                    ("Therapeutic", "antibodyName"),
                    ("Format", "antiBodyFormat"),
                    ("CH1 Isotype", "ch1Isotype"),
                    ("VD LC", "VD_LC"),
                    (clinTrialCol, "maxClinicalPhase"),
                    ("Est. Status", "status"),
                    ("Target", "target"),
                    ("Conditions Approved", "conditionsApproved"),
                    ("Conditions Active", "conditionsActive"),
                ]:
                    if kTup[0] in rD and rD[kTup[0]] not in ["na", "na;na"]:
                        qD[kTup[1]] = rD[kTup[0]]
                    else:
                        qD[kTup[1]] = None
                    if kTup[0] not in rD:
                        logger.error("SabDab key %r missing in input dataset %r", kTup[0], list(rD.keys()))
                tD[rD["Therapeutic"]] = qD
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
        """Fetch and read SAbDab2 antibody assignment data (JSON from the rcsb-pdb-annotations endpoint).

        The source contains one record per heavy/light chain pairing, with nested chain and antigen details.
        Records are collapsed to one assignment per antibody chain, since a chain can appear in several pairings.

        Args:
            dirPath (str): local directory in which to store the downloaded data

        Returns:
            dict: assignments keyed by '<entryId>.<authAsymId>' (entry ID lower case)

        Raises:
            ValueError: if the data cannot be fetched or yields no assignments
        """
        startTime = time.time()
        targetUrl = kwargs.get("assignmentUrl", "https://sabdab.opig.stats.ox.ac.uk/api/rcsb-pdb-annotations")
        fU = FileUtil()
        fU.mkdir(dirPath)
        dumpPath = os.path.join(dirPath, "sabdab_summary_all.json")
        logger.info("Fetching url %s path %s", targetUrl, dumpPath)
        ok = fU.get(targetUrl, dumpPath)
        if not ok:
            raise ValueError("Fetching failed for SAbDab assignment data")
        rDL = self.__mU.doImport(dumpPath, fmt="json")
        logger.info("SAbDab raw records (%d)", len(rDL))
        #
        # NOTE: The following keys from the old SAbDab summary file are no longer provided (and not populated here):
        #   Hchain/Lchain (the auth_asym_id is now used directly as part of the assignment key),
        #   antigen_chain, antigen_type, antigen_het_name.
        # antigen_het_name was dropped by SAbDab2 on purpose (inconsistent meaning); antigen_chain and antigen_type
        # have no direct equivalent (antigen instances only carry "entity type" and chain IDs of the antigen itself).
        # None of these were used by downstream consumers (e.g., dictionary utils ETL).
        #
        # Per-chain accumulators, so chains seen in multiple pairings are merged rather than overwritten
        chainD = {}
        for rD in rDL:
            entryId = (rD.get("PDB ID") or "").lower()
            if not entryId:
                continue
            # NOTE: All antigen instances in the source data are retained as-is, including ions and haptens
            # (e.g., "CHLORIDE ION", "GLYCEROL", "1,2-ETHANEDIOL"), even where they are crystallization additives.
            antigenNameL = [antD["antigen name"] for antD in rD.get("antigen_instances") or [] if antD.get("antigen name")]
            for chainKey in ["heavy chain", "light chain"]:
                cD = rD.get(chainKey)
                authAsymId = cD.get("PDB auth_asym_id") if cD else None
                if not authAsymId:
                    continue
                asgD = chainD.setdefault(entryId + "." + authAsymId, {"pdb": entryId, "model": rD.get("model"), "antigens": set()})
                asgD["antigens"].update(antigenNameL)
                # V gene subgroup, e.g. "IGLV1 (Homsap),IGKV3 (Homsap)" -> "IGLV1,IGKV3"; many chains have none (IMGT-derived)
                subgroup = self.__cleanSubgroup(cD.get("V gene subgroup"))
                if subgroup:
                    asgD.setdefault("heavy_subclass" if chainKey == "heavy chain" else "light_subclass", subgroup)
                if chainKey == "light chain":
                    cType = {"λ": "Lambda", "κ": "Kappa"}.get(cD.get("chain type"), cD.get("chain type"))
                    if cType:
                        asgD.setdefault("light_ctype", cType)
        #
        aD = {}
        for instId, cD in chainD.items():
            antigens = cD.pop("antigens")
            if antigens:
                cD["antigen_name"] = " | ".join(sorted(antigens))
            aD[instId] = cD
            # Also index extended entry IDs (pdb_0000xxxx) by their short ID, for lookups that only have the short ID
            if instId.startswith("pdb_0000"):
                aD[instId.removeprefix("pdb_0000")] = cD
        #
        logger.info("Fetched (%d) SAbDab assignment records.", len(aD))
        if not aD:
            raise ValueError("No SAbDab assignments found in %r" % dumpPath)
        logger.info("Completed reload at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return aD

    @staticmethod
    def __cleanSubgroup(value):
        """Strip parenthetical organism names and duplicates from a V gene subgroup string."""
        if not value:
            return None
        gL = []
        for g in re.sub(r"\s*\([^)]*\)", "", value).split(","):
            g = g.strip()
            if g and g not in gL:
                gL.append(g)
        return ",".join(gL) or None

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
                hsKey = "Heavy Sequence" if "Heavy Sequence" in rD else "HeavySequence"
                hSeq = rD[hsKey] if rD[hsKey] != "na" else None
                if not hSeq:
                    continue
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
                lsKey = "Light Sequence" if "Light Sequence" in rD else "LightSequence"
                lSeq = rD[lsKey] if rD[lsKey] != "na" else None
                if not lSeq:
                    continue
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
