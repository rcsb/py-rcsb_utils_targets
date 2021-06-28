##
#  File:           PharosTargetProvider.py
#  Date:           11-Jun-2021 jdw
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
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.seq.UniProtIdMappingProvider import UniProtIdMappingProvider

logger = logging.getLogger(__name__)


class PharosTargetProvider(StashableBase):
    """Accessors for Pharos target assignments."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirName = "Pharos-targets"
        super(PharosTargetProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        reloadDb = kwargs.get("reloadDb", False)
        fromDb = kwargs.get("fromDb", False)
        useCache = kwargs.get("useCache", False)
        pharosDumpUrl = kwargs.get("pharosDumpUrl", None)
        mysqlUser = kwargs.get("mysqlUser", None)
        mysqlPassword = kwargs.get("mysqlPassword", None)
        self.__version = None
        if reloadDb or fromDb:
            self.__reload(self.__dirPath, reloadDb=reloadDb, fromDb=fromDb, useCache=useCache, pharosDumpUrl=pharosDumpUrl, mysqlUser=mysqlUser, mysqlPassword=mysqlPassword)
        #

    def testCache(self):
        return True

    def getVersion(self):
        return self.__version

    def __reload(self, dirPath, reloadDb=False, fromDb=False, useCache=False, pharosDumpUrl=None, mysqlUser=None, mysqlPassword=None):
        startTime = time.time()
        pharosSelectedTables = ["drug_activity", "cmpd_activity", "target", "protein", "t2tc"]
        pharosDumpUrl = pharosDumpUrl if pharosDumpUrl else "http://juniper.health.unm.edu/tcrd/download/latest.sql.gz"
        pharosReadmeUrl = "http://juniper.health.unm.edu/tcrd/download/latest.README"
        ok = False
        fU = FileUtil()
        pharosDumpFileName = fU.getFileName(pharosDumpUrl)
        pharosDumpPath = os.path.join(dirPath, pharosDumpFileName)
        pharosUpdatePath = os.path.join(dirPath, "pharos-update.sql")
        pharosReadmePath = os.path.join(dirPath, "pharos-readme.txt")
        logPath = os.path.join(dirPath, "pharosLoad.log")
        #
        fU.mkdir(dirPath)
        #

        exU = ExecUtils()
        #
        if reloadDb:
            logger.info("useCache %r pharosDumpPath %r", useCache, pharosDumpPath)
            if useCache and self.__mU.exists(pharosDumpPath):
                ok = True
            else:
                logger.info("Fetching url %s path %s", pharosDumpUrl, pharosDumpPath)
                ok1 = fU.get(pharosDumpUrl, pharosDumpPath)
                ok2 = fU.get(pharosReadmeUrl, pharosReadmePath)
                logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok1 and ok2, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            # ---
            readmeLines = self.__mU.doImport(pharosReadmePath, fmt="list")
            self.__version = readmeLines[0].split(" ")[1][1:] if readmeLines else "6"
            # ---
            logger.info("Filtering SQL dump %r for selected tables %r", pharosDumpFileName, pharosSelectedTables)
            doWrite = True
            # Note: the pharos dump file latest.sql.gz is not gzipped
            with open(pharosDumpPath, "r") as ifh, open(pharosUpdatePath, "w") as ofh:
                for line in ifh:
                    if line.startswith("-- Table structure for table"):
                        tN = line.split(" ")[-1][1:-2]
                        doWrite = True if tN in pharosSelectedTables else False
                    if doWrite:
                        ofh.write(line)
            # ---
            ok = exU.run(
                "mysql",
                execArgList=["-v", "-u", mysqlUser, "--password=%s" % mysqlPassword, "-e", "create database if not exists tcrd6;"],
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
            shellCmd = 'trap "" SIGHUP SIGINT SIGTERM; nohup mysql -u %s --password=%s tcrd6 < %s >& %s' % (mysqlUser, mysqlPassword, pharosUpdatePath, logPath)
            ok = exU.runShell(
                shellCmd,
                outPath=None,
                inpPath=None,
                outAppend=True,
                timeOut=None,
            )
            logger.info("SQL dump restore status %r", ok)
        # --
        if fromDb:
            for tbl in pharosSelectedTables:
                outPath = os.path.join(dirPath, "%s.tdd" % tbl)
                # if useCache and self.__mU.exists(outPath):
                #   continue
                ok = exU.run(
                    "mysql",
                    execArgList=["-u", mysqlUser, "--password=%s" % mysqlPassword, "-e", "use tcrd6; select * from %s;" % tbl],
                    outPath=outPath,
                    outAppend=False,
                    timeOut=None,
                    suppressStderr=True,
                )
                logger.info("SQL table %s export status %r", tbl, ok)
        return ok

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
