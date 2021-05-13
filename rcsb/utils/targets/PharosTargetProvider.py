##
#  File:           PharosTargetProvider.py
#  Date:           9-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for Pharos target assignments.

"""

import logging
import os.path
import time

from rcsb.utils.io.ExecUtils import ExecUtils
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashUtil import StashUtil
from rcsb.utils.seq.UniProtIdMappingProvider import UniProtIdMappingProvider

logger = logging.getLogger(__name__)


class PharosTargetProvider:
    """Accessors for Pharos target assignments."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirName = "Pharos-targets"
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        reloadDb = kwargs.get("reloadDb", False)
        fromDb = kwargs.get("fromDb", False)
        if reloadDb or fromDb:
            self.__reload(self.__dirPath, reloadDb=reloadDb, fromDb=fromDb, **kwargs)
        #

    def testCache(self):
        return True

    def __reload(self, dirPath, reloadDb=False, fromDb=False, **kwargs):
        startTime = time.time()
        useCache = kwargs.get("useCache", True)
        pharosDumpUrl = kwargs.get("pharosDumpUrl", "http://juniper.health.unm.edu/tcrd/download/latest.sql.gz")

        ok = False
        fU = FileUtil()
        pharosDumpFileName = fU.getFileName(pharosDumpUrl)
        pharosDumpPath = os.path.join(dirPath, pharosDumpFileName)
        logPath = os.path.join(dirPath, "pharosLoad.log")
        #
        fU.mkdir(dirPath)
        #
        mysqlUser = kwargs.get("mysqlUser", None)
        mysqlPassword = kwargs.get("mysqlPassword", None)
        exU = ExecUtils()
        #
        if reloadDb:
            logger.info("useCache %r pharosDumpPath %r", useCache, pharosDumpPath)
            if useCache and self.__mU.exists(pharosDumpPath):
                ok = True
            else:
                logger.info("Fetching url %s path %s", pharosDumpUrl, pharosDumpPath)
                ok = fU.get(pharosDumpUrl, pharosDumpPath)
                logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            # ---
            ok = exU.run(
                "mysql",
                execArgList=["-v", "-u", mysqlUser, "--password=%s" % mysqlPassword, "-e", "drop database if exists tcrd6; create database tcrd6;"],
                outPath=logPath,
                outAppend=False,
                timeOut=None,
            )
            # ok = exU.run(
            #     "mysql",
            #     execArgList=["-u", mysqlUser, "--password=%s" % mysqlPassword, "tcrd6"],
            #     outPath=logPath,
            #     inpPath=pharosDumpPath,
            #     outAppend=True,
            #     timeOut=None,
            # )
            shellCmd = 'trap "" SIGHUP SIGINT SIGTERM; nohup mysql -u %s --password=%s tcrd6 < %s >& %s' % (mysqlUser, mysqlPassword, pharosDumpPath, logPath)
            ok = exU.runShell(
                shellCmd,
                outPath=None,
                inpPath=None,
                outAppend=True,
                timeOut=None,
            )
        # --
        if fromDb:
            for tbl in ["drug_activity", "cmpd_activity", "target", "protein", "t2tc"]:
                outPath = os.path.join(dirPath, "%s.tdd" % tbl)
                if useCache and self.__mU.exists(outPath):
                    continue
                ok = exU.run(
                    "mysql",
                    execArgList=["-u", mysqlUser, "--password=%s" % mysqlPassword, "-e", "use tcrd6; select * from %s;" % tbl],
                    outPath=os.path.join(dirPath, "%s.tdd" % tbl),
                    outAppend=False,
                    timeOut=None,
                    suppressStderr=True,
                )
        return True

    def exportProteinFasta(self, fastaPath, taxonPath, addTaxonomy=False):
        try:
            proteinFilePath = os.path.join(self.__dirPath, "protein.tdd")
            pDL = self.__mU.doImport(proteinFilePath, fmt="tdd", rowFormat="dict")
            fD = {}
            taxonL = []
            if addTaxonomy:
                umP = UniProtIdMappingProvider(self.__cachePath)
                umP.reload(useCache=True)
                #
                for pD in pDL:
                    unpId = pD["uniprot"]
                    proteinId = pD["id"]
                    seq = pD["seq"]
                    taxId = umP.getMappedId(unpId, mapName="NCBI-taxon")
                    taxId = taxId if taxId else "-1"
                    cD = {"sequence": seq, "uniprotId": unpId, "proteinId": proteinId, "taxId": taxId}
                    seqId = ""
                    cL = []
                    for k, v in cD.items():
                        if k in ["sequence"]:
                            continue
                        cL.append(str(v))
                        cL.append(str(k))
                    seqId = "|".join(cL)
                    fD[seqId] = cD
                    taxonL.append("%s\t%s" % (seqId, taxId))
                ok = self.__mU.doExport(taxonPath, taxonL, fmt="list")
            else:
                for pD in pDL:
                    unpId = pD["uniprot"]
                    proteinId = pD["id"]
                    seq = pD["seq"]
                    cD = {"sequence": seq, "uniprotId": unpId, "proteinId": proteinId}
                    seqId = ""
                    cL = []
                    for k, v in cD.items():
                        if k in ["sequence"]:
                            continue
                        cL.append(str(v))
                        cL.append(str(k))
                    seqId = "|".join(cL)
                    fD[seqId] = cD
            #
            logger.info("Writing %d pharos targets to %s", len(fD), fastaPath)
            ok = self.__mU.doExport(fastaPath, fD, fmt="fasta", makeComment=True)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def exportCofactors(self, cofactorDataPath, fmt="json"):
        targetD = {}
        cofactorFilePath = os.path.join(self.__dirPath, "drug_activity.tdd")
        cfDL = self.__mU.doImport(cofactorFilePath, fmt="tdd", rowFormat="dict")
        targetD = self.__extactCofactorData(cfDL)
        #
        cofactorFilePath = os.path.join(self.__dirPath, "cmpd_activity.tdd")
        cfDL = self.__mU.doImport(cofactorFilePath, fmt="tdd", rowFormat="dict")
        targetD.update(self.__extactCofactorData(cfDL))
        #
        tD = self.__getTargetDetails()
        for tId in targetD:
            if tId in tD:
                for k, v in tD[tId].items():
                    targetD[tId][k] = v
        #
        ok = self.__mU.doExport(cofactorDataPath, targetD, fmt=fmt)
        return ok

    def __getTargetDetails(self):
        # protein.tdd
        # id	name	description	uniprot	up_version	geneid	sym	family	chr	seq	dtoid	stringid	dtoclass
        rD = {}
        try:
            proteinFilePath = os.path.join(self.__dirPath, "protein.tdd")
            pDL = self.__mU.doImport(proteinFilePath, fmt="tdd", rowFormat="dict")
            for pD in pDL:
                proteinId = pD["id"]
                unpId = pD["uniprot"] if "uniprot" in pD and pD["uniprot"] != "NULL" else None
                descr = pD["description"] if "description" in pD and pD["description"] != "NULL" else None
                geneId = pD["geneid"] if "geneid" in pD and pD["geneid"] != "NULL" else None
                dtoId = pD["dtoid"] if "dtoid" in pD and pD["dtoid"] != "NULL" else None
                dtoClass = pD["dtoclass"] if "dtoclass" in pD and pD["dtoclass"] != "NULL" else None
                rD[proteinId] = {"unpId": unpId, "name": descr, "geneId": geneId, "dtoId": dtoId, "dtoClass": dtoClass}
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return rD

    def __extactCofactorData(self, cfDL):
        """Extra ids, activity and moa data for drugs and cofactors.

        Args:
            cfDL (list): list of dictionaries of containing pharos exported db data.

        Returns:
            dict: dictionary of extract cofactor data
        """
        try:
            qD = {}
            targetD = {}
            for cfD in cfDL:
                tId = cfD["target_id"]
                qD = {}
                qD["smiles"] = cfD["smiles"] if "smiles" in cfD and cfD["smiles"] != "N" else None
                qD["chemblId"] = cfD["cmpd_chemblid"] if "cmpd_chemblid" in cfD else None
                qD["chemblId"] = cfD["cmpd_id_in_src"] if "catype" in cfD and cfD["catype"].upper() == "CHEMBL" else qD["chemblId"]
                qD["pubChemId"] = cfD["cmpd_pubchem_cid"] if "cmpd_pubchem_cid" in cfD else None
                #
                qD["activity"] = cfD["act_value"] if "act_value" in cfD and cfD["act_value"] != "NULL" else None
                qD["activityType"] = cfD["act_type"] if "act_type" in cfD and cfD["act_type"] != "NULL" else None
                if qD["activity"] is not None:
                    qD["activity"] = float(qD["activity"])
                #
                qD["moa"] = cfD["action"] if "action" in cfD and cfD["moa"] == "1" else None
                qD["name"] = cfD["drug"] if "drug" in cfD else None
                qD["name"] = cfD["cmpd_name_in_src"] if "cmpd_name_in_src" in cfD and cfD["cmpd_name_in_src"] != "NULL" else qD["name"]
                targetD[tId] = {ky: qD[ky] for ky in qD if qD[ky] is not None}
            #
        except Exception as e:
            logger.exception("Failing with %r %s", qD, str(e))
        return targetD

    def restore(self, cfgOb, configName):
        ok = False
        try:
            startTime = time.time()
            url = cfgOb.get("STASH_SERVER_URL", sectionName=configName)
            userName = cfgOb.get("_STASH_AUTH_USERNAME", sectionName=configName)
            password = cfgOb.get("_STASH_AUTH_PASSWORD", sectionName=configName)
            basePath = cfgOb.get("_STASH_SERVER_BASE_PATH", sectionName=configName)
            ok = self.__fromStash(url, basePath, userName=userName, password=password)
            logger.info("Recovered Pharos data file from stash (%r)", ok)
            if not ok:
                urlFallBack = cfgOb.get("STASH_SERVER_FALLBACK_URL", sectionName=configName)
                ok = self.__fromStash(urlFallBack, basePath, userName=userName, password=password)
                logger.info("Recovered Pharos data file from fallback stash (%r)", ok)
            #
            logger.info("Completed recovery (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return ok

    def backup(self, cfgOb, configName):
        ok1 = ok2 = False
        try:
            startTime = time.time()
            userName = cfgOb.get("_STASH_AUTH_USERNAME", sectionName=configName)
            password = cfgOb.get("_STASH_AUTH_PASSWORD", sectionName=configName)
            basePath = cfgOb.get("_STASH_SERVER_BASE_PATH", sectionName=configName)
            url = cfgOb.get("STASH_SERVER_URL", sectionName=configName)
            urlFallBack = cfgOb.get("STASH_SERVER_FALLBACK_URL", sectionName=configName)
            ok1 = self.__toStash(url, basePath, userName=userName, password=password)
            ok2 = self.__toStash(urlFallBack, basePath, userName=userName, password=password)
            logger.info("Completed Pharos backup (%r/%r) at %s (%.4f seconds)", ok1, ok2, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok1 & ok2

    def __toStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        """Copy tar and gzipped bundled cache data to remote server/location.

        Args:
            url (str): server URL (e.g. sftp://hostname.domain) None for local host
            stashRemoteDirPath (str): path to target directory on remote server
            userName (str, optional): server username. Defaults to None.
            password (str, optional): server password. Defaults to None.
            remoteStashPrefix (str, optional): channel prefix. Defaults to None.

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            stU = StashUtil(os.path.join(self.__cachePath, "stash"), self.__dirName)
            ok = stU.makeBundle(self.__cachePath, [self.__dirName])
            if ok:
                ok = stU.storeBundle(url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.error("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok

    def __fromStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        """Restore local cache from a tar and gzipped bundle to fetched from a remote server/location.

        Args:
            url (str): server URL (e.g. sftp://hostname.domain) None for local host
            stashRemoteDirPath (str): path to target directory on remote server
            userName (str, optional): server username. Defaults to None.
            password (str, optional): server password. Defaults to None.
            remoteStashPrefix (str, optional): channel prefix. Defaults to None.

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            stU = StashUtil(os.path.join(self.__cachePath, "stash"), self.__dirName)
            ok = stU.fetchBundle(self.__cachePath, url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.error("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok
