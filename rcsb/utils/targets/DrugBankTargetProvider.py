##
# File: DrugBankProvider.py
# Author:  J. Westbrook
# Date:  8-Nov-2020
#
# Manage access to DrugBank target data.
#
# Update:
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

# from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class DrugBankTargetProvider(object):
    """Utilities to manage DrugBank target FASTA data."""

    def __init__(self, **kwargs):
        #
        self.__fastaPathList = self.__reload(**kwargs)
        self.__version = None

    def testCache(self):
        try:
            return True
        except Exception as e:
            logger.debug("Failing with %s", str(e))
        return False

    def getVersion(self):
        return None

    def __reload(self, **kwargs):
        """Reload DrugBank target FASTA data files.

        Args:
            cachePath (str, optional): path to top cache directory
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

        dirPath = os.path.join(kwargs.get("cachePath", "."), "DrugBank-targets")
        useCache = kwargs.get("useCache", True)
        username = kwargs.get("username", None)
        password = kwargs.get("password", None)
        # mU = MarshalUtil(workPath=dirPath)
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
        #
        return retFilePathList
