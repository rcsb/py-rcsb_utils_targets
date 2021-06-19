##
# File: DrugBankProvider.py
# Author:  J. Westbrook
# Date:  8-Nov-2020
#
# Manage access to DrugBank target data.
#
# Update:
#
#
##

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import datetime
import logging
import os
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seq.UniProtIdMappingProvider import UniProtIdMappingProvider

logger = logging.getLogger(__name__)


class DrugBankTargetProvider(object):
    """Utilities to manage DrugBank target FASTA data.

    UniProt mapping URL template:

        https://go.drugbank.com/polypeptides/P23458
    """

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "DrugBank-targets")
        self.__version = None
        self.__cfD = None
        self.__fastaPathList = self.__reloadFasta(self.__dirPath, **kwargs)
        self.__mU = MarshalUtil(workPath=self.__dirPath)

    def testCache(self):
        try:
            return len(self.__fastaPathList) >= 4
        except Exception as e:
            logger.debug("Failing with %s", str(e))
        return False

    def getAssignmentVersion(self):
        return self.__version if self.__version else datetime.datetime.now().strftime("%Y-%m-%d")

    def __reloadFasta(self, dirPath, **kwargs):
        """Reload DrugBank target FASTA data files.

        Args:
            dirPath (str, optional): path to DrugBank cache directory
            useCache (bool, optional): flag to use cached files. Defaults to True.

        Returns:

        """
        startTime = time.time()
        logger.info("Starting db reload at %s", time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        retFilePathList = []
        urlTargetL = [
            "https://go.drugbank.com/releases/latest/downloads/target-all-polypeptide-sequences",
            "https://go.drugbank.com/releases/latest/downloads/enzyme-all-polypeptide-sequences",
            "https://go.drugbank.com/releases/latest/downloads/carrier-all-polypeptide-sequences",
            "https://go.drugbank.com/releases/latest/downloads/transporter-all-polypeptide-sequences",
        ]

        useCache = kwargs.get("useCache", True)
        username = kwargs.get("username", None)
        password = kwargs.get("password", None)
        #
        if not username or not password:
            return retFilePathList
        #
        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        if not useCache:
            #  Clear any cached files
            for urlTarget in urlTargetL:
                baseFileName = fU.getFileName(urlTarget)
                zipFileName = baseFileName + ".fasta.zip"
                retFileName = baseFileName + ".fa"
                for fn in [baseFileName, zipFileName, retFileName]:
                    try:
                        fp = os.path.join(dirPath, fn)
                        os.remove(fp)
                    except Exception:
                        pass
        #
        ok = False
        if useCache:
            ok = True
            for urlTarget in urlTargetL:
                baseFileName = fU.getFileName(urlTarget)
                retFileName = baseFileName + ".fa"
                retFilePath = os.path.join(dirPath, retFileName)
                ok = fU.exists(retFilePath)
                if not ok:
                    break
                retFilePathList.append(retFilePath)
        #
        logger.info("Using cached files %r", ok)
        if not useCache or not ok:
            if not username or not password:
                logger.warning("Missing credentials for DrugBank file download...")

            for urlTarget in urlTargetL:
                baseFileName = fU.getFileName(urlTarget)
                zipFileName = baseFileName + ".fasta.zip"
                retFileName = baseFileName + ".fa"
                zipFilePath = os.path.join(dirPath, zipFileName)
                retFilePath = os.path.join(dirPath, retFileName)
                basePath = os.path.join(dirPath, baseFileName)
                logger.info("Fetching url %s for FASTA target file %s", urlTarget, baseFileName)
                ok = fU.get(urlTarget, zipFilePath, username=username, password=password)
                endTime = time.time()
                logger.info("Completed db fetch at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
                #
                ok = fU.unbundleZipfile(zipFilePath, dirPath=basePath)
                fU.put(os.path.join(basePath, "protein.fasta"), retFilePath)
                endTime = time.time()
                logger.info("Completed unzip at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
                retFilePathList.append(retFilePath)
        return retFilePathList

    def exportFasta(self, fastaPath, taxonPath, addTaxonomy=False):
        ok = self.__consolidateFasta(fastaPath, taxonPath, self.__fastaPathList, addTaxonomy=addTaxonomy)
        return ok

    def __consolidateFasta(self, fastaPath, taxonPath, inpPathList, addTaxonomy=False):
        #
        drugBankTargetMapPath = self.__getTargetDrugMapPath()
        #
        oD = {}
        uD = {}
        taxonL = []
        try:
            if addTaxonomy:
                umP = UniProtIdMappingProvider(self.__cachePath)
                umP.reload(useCache=True)
            #
            for fp in inpPathList:
                fD = self.__mU.doImport(fp, fmt="fasta", commentStyle="default")
                for seqId, sD in fD.items():
                    tL = seqId[seqId.find("(") + 1 : seqId.find(")")]
                    dbIdL = [v.strip() for v in tL.split(";")]
                    seq = sD["sequence"]
                    tL = seqId.split("|")
                    unpId = tL[1].split(" ")[0]
                    cD = {"sequence": seq, "uniprotId": unpId}
                    if addTaxonomy:
                        taxId = umP.getMappedId(unpId, mapName="NCBI-taxon")
                        cD["taxId"] = taxId if taxId else -1
                    #
                    seqId = ""
                    cL = []
                    for k, v in cD.items():
                        if k in ["sequence"]:
                            continue
                        cL.append(str(v))
                        cL.append(str(k))
                    seqId = "|".join(cL)
                    oD[seqId] = cD
                    if addTaxonomy:
                        taxonL.append("%s\t%s" % (seqId, taxId))
                    uD.setdefault(unpId, []).extend(dbIdL)

            ok1 = self.__mU.doExport(fastaPath, oD, fmt="fasta", makeComment=True)
            tS = datetime.datetime.now().isoformat()
            vS = datetime.datetime.now().strftime("%Y-%m-%d")
            ok2 = self.__mU.doExport(drugBankTargetMapPath, {"version": vS, "created": tS, "cofactors": uD}, fmt="json")
            ok3 = True
            if addTaxonomy:
                ok3 = self.__mU.doExport(taxonPath, taxonL, fmt="list")
            return ok1 & ok2 & ok3
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return False

    def __getTargetDrugMapPath(self):
        return os.path.join(self.__dirPath, "drugbank_target_drug_map.json")

    def hasCofactor(self, unpId):
        if not self.__cfD:
            self.__cfD = self.__reloadCofactors()
        try:
            return unpId in self.__cfD
        except Exception:
            return False

    def getCofactors(self, unpId):
        if not self.__cfD:
            self.__cfD = self.__reloadCofactors()
        try:
            return self.__cfD[unpId]
        except Exception:
            return []

    def __reloadCofactors(self):
        try:
            qD = self.__mU.doImport(self.__getTargetDrugMapPath(), fmt="json")
            self.__version = qD["version"]
            return qD["cofactors"]
        except Exception as e:
            logger.error("Failing with %s", str(e))
        return {}
