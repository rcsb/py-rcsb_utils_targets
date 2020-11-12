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

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"


import logging
import os
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seqalign.MiscUtils import MiscUtils

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
        self.__fastaPathList = self.__reload(self.__dirPath, **kwargs)
        self.__consolidateFasta(self.__dirPath, self.__fastaPathList, **kwargs)
        self.__version = None

    def testCache(self):
        try:
            return True
        except Exception as e:
            logger.debug("Failing with %s", str(e))
        return False

    def getVersion(self):
        return None

    def __reload(self, dirPath, **kwargs):
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

    def __consolidateFasta(self, dirPath, inpPathList, **kwargs):
        #
        drugBankfastaFileName = kwargs.get("drugBankfastaFileName ", "drugbank_targets.fa")
        drugBankTargetMapFileName = kwargs.get("drugBankTargetMapFileName", "drugbank_target_drug_map.json")
        addTaxonomy = kwargs.get("addTaxonomy", False)
        #
        mU = MarshalUtil(workPath=dirPath)
        oD = {}
        uD = {}
        try:
            if addTaxonomy:
                miscU = MiscUtils()
                outDirPath = os.path.join(self.__cachePath, "uniprot_id_mapping_selected")
                taxMapFileName = "uniprot_taxonomy.pic"
                taxD = miscU.getUniprotXref(13, outDirPath, taxMapFileName, fmt="pickle", useCache=True)
            #
            for fp in inpPathList:
                fD = mU.doImport(fp, fmt="fasta", commentStyle="default")
                for seqId, sD in fD.items():
                    tL = seqId[seqId.find("(") + 1 : seqId.find(")")]
                    dbIdL = [v.strip() for v in tL.split(";")]
                    seq = sD["sequence"]
                    tL = seqId.split("|")
                    unpId = tL[1].split(" ")[0]
                    if addTaxonomy:
                        taxId = taxD[unpId] if unpId in taxD else "-1"
                        seqId = unpId + "|" + taxId
                    else:
                        seqId = unpId
                    uD.setdefault(unpId, []).extend(dbIdL)
                    oD[seqId] = {"sequence": seq}
            ok1 = mU.doExport(os.path.join(dirPath, drugBankfastaFileName), oD, fmt="fasta")
            ok2 = mU.doExport(os.path.join(dirPath, drugBankTargetMapFileName), uD, fmt="json")
            return ok1 & ok2
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return False
