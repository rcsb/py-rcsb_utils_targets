##
#  File:           SAbDabTargetProvider.py
#  Date:           18-Jun-2021 jdw
#
#  Updated:
#   9-Feb-2023 aae  Find Highest_Clin_Trial column regardless of month
#   1-Jul-2024 dwp  Update SAbDab data parsing (following change in source data headers)
#  22-Jul-2024 dwp  Re-update SAbDab data parsing (following change in source data headers)
#   2-Jul-2026 dwp  Use backup SAbDab summary file and TheraSAbDab file during transition phase to SAbDab2
##
"""
Accessors for Thera-SAbDab(Therapeutic Structural Antibody Database) target data.
"""

import datetime
import logging
import os.path
import time
import json
import re
import requests

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
            featureKey (str): assignment feature key: one of pdb|Hchain|Lchain|model|antigen_chain|antigen_type|
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
        targetUrl = kwargs.get(
            "targetUrl",
            "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv"
        )
        targetFallbackUrl = "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets_stash/refs/heads/development/stash/SAbDab-backup/TheraSAbDab_SeqStruc_OnlineDownload.csv"
        #
        # TODO: need to update UI links to SabDab
        # e.g., https://www.rcsb.org/annotations/1BEY#antibodyTheraSAbDab
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
                logger.error("Fetching failed for Thera-SAbDab target data. Trying fallback.")
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
        """Fetch and read SAbDab antibody assignment data.

        The SAbDab JSON format contains one record per PDB/model and
        nested heavy/light chain and antigen instance information.

        Args:
            dirPath (str): Local directory in which to store the downloaded data.

        Returns:
            dict: Assignment records keyed by '<pdb_id>.<auth_asym_id>'.
        """
        startTime = time.time()
        aD = {}
        ok = False

        try:
            targetUrl = kwargs.get("assignmentUrl", "https://sabdab.opig.stats.ox.ac.uk/api/rcsb-pdb-annotations")
            # TODO: create a fallback JSON file and post to stash

            fU = FileUtil()
            dumpFileName = "sabdab_summary_all.json"
            fU.mkdir(dirPath)
            dumpPath = os.path.join(dirPath, dumpFileName)
            logger.info("Fetching url %s path %s", targetUrl, dumpPath)
            response = requests.get(targetUrl, headers={"Accept-Encoding": "gzip"}, timeout=60)
            response.raise_for_status()
            rDL = response.json()
            with open(dumpPath, "w", encoding="utf-8") as f:
                json.dump(rDL, f, indent=2, ensure_ascii=False)
            logger.info("SAbDab raw records (%d)", len(rDL))
            if rDL:
                logger.debug("rD keys %r", list(rDL[0].keys()))

            for rD in rDL:
                # Convert extended ID "pdb_00009m5h" -> short ID "9m5h"
                # TODO: Stop doing this conversion when switching to Beta-Archive loading for extended IDs
                pdbId = rD.get("PDB ID")
                if pdbId and pdbId.startswith("pdb_0000"):
                    pdbId = pdbId.removeprefix("pdb_0000")
                if not pdbId:
                    continue

                model = rD.get("model")

                # Collect antigen information
                antigenInstances = rD.get("antigen_instances") or []
                antigenNames = []

                for antigenD in antigenInstances:
                    antigenName = antigenD.get("antigen name")
                    if antigenName:
                        antigenNames.append(str(antigenName))

                # Remove duplicate values while preserving their order.
                antigenNames.sort()
                antigenNames = list(set(antigenNames))
                antigenName = " | ".join(antigenNames) if antigenNames else None

                commonD = {
                    "pdb": pdbId,
                    "model": model,
                    "antigen_name": antigenName,
                }

                # Heavy-chain assignment
                heavyD = rD.get("heavy chain")
                if heavyD:
                    authAsymIdH = heavyD.get("PDB auth_asym_id")
                    if authAsymIdH:
                        assignmentD = dict(commonD)
                        # TODO: Change labeling of "subclass" to "subgroup" on Annotations UI page
                        # Unfortunately, can't easily change the internal labeling since it's defined in the schema ("SABDAB_ANTIBODY_LIGHT_CHAIN_SUBCLASS")
                        heavySubclass = heavyD.get("V gene subgroup")
                        if heavySubclass:
                            # remove parenthetical organisms (e.g., "IGLV1 (Homsap),IGKV3 (Homsap)")
                            heavySubclassExtract = re.sub(r"\s*\([^)]*\)", "", heavySubclass)
                            assignmentD["heavy_subclass"] = heavySubclassExtract

                        aD[pdbId + "." + authAsymIdH] = assignmentD

                # Light-chain assignment
                lightD = rD.get("light chain")
                if lightD:
                    authAsymIdL = lightD.get("PDB auth_asym_id")
                    if authAsymIdL:
                        assignmentD = dict(commonD)
                        lightSubclass = lightD.get("V gene subgroup")
                        if lightSubclass:
                            # remove parenthetical organisms (e.g., "IGLV1 (Homsap),IGKV3 (Homsap)")
                            lightSubclassExtract = re.sub(r"\s*\([^)]*\)", "", lightSubclass)
                            assignmentD["light_subclass"] = lightSubclassExtract

                        lightCtype = lightD.get("chain type")
                        if lightCtype:
                            if lightCtype == "λ":
                                lightCtype = "Lambda"
                            elif lightCtype == "κ":
                                lightCtype = "Kappa"
                            assignmentD["light_ctype"] = lightCtype

                        aD[pdbId + "." + authAsymIdL] = assignmentD

            logger.info("Fetched (%d) SAbDab assignment records.", len(aD))
            ok = len(aD) > 0

        except Exception as e:
            logger.exception("Failing with %s", str(e))

        logger.info(
            "Completed reload (%r) at %s (%.4f seconds)",
            ok,
            time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
            time.time() - startTime,
        )

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
