##
# File: DrugBankProvider.py
# Author:  J. Westbrook
# Date:  8-Nov-2020
#
# Manage access to DrugBank target data.
#
# Update:
#  22-Jun-2021 jdw replaced the download of FASTA files with processing the main DrugBank download data file.
#                  Leaving the deprecated FASTA code for now -
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

from rcsb.utils.chemref.DrugBankProvider import DrugBankProvider
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seq.UniProtIdMappingProvider import UniProtIdMappingProvider

logger = logging.getLogger(__name__)


class DrugBankTargetProvider(object):
    """Utilities to manage DrugBank target FASTA data.

    Note: UniProt mapping URL template:

        https://go.drugbank.com/polypeptides/P23458
    """

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "DrugBank-targets")
        self.__cfD = None
        # self.__fastaPathList = self.__reloadFasta(self.__dirPath, **kwargs)
        self.__dbP = self.__reloadDrugBankProvider(**kwargs)
        self.__version = self.__dbP.getVersion() if self.__dbP else None
        self.__mU = MarshalUtil(workPath=self.__dirPath)

    def testCache(self):
        try:
            # return len(self.__fastaPathList) >= 4
            return self.__dbP is not None
        except Exception as e:
            logger.debug("Failing with %s", str(e))
        return False

    def getAssignmentVersion(self):
        return self.__version if self.__version else datetime.datetime.now().strftime("%Y-%m-%d")

    def getFastaPath(self):
        return os.path.join(self.__cachePath, "FASTA", "drugbank")

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
        # ok = self.__consolidateFasta(fastaPath, taxonPath, self.__fastaPathList, addTaxonomy=addTaxonomy)
        ok = self.__buildResourceFiles(fastaPath, taxonPath, addTaxonomy=addTaxonomy)
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

    def __reloadDrugBankProvider(self, **kwargs):
        useCache = kwargs.get("useCache", True)
        un = kwargs.get("username", None)
        pw = kwargs.get("password", None)
        return DrugBankProvider(cachePath=self.__cachePath, useCache=useCache, username=un, password=pw)

    def __buildResourceFiles(self, fastaPath, taxonPath, addTaxonomy=False):
        """Build DrugBank FASTA and resource files from the full DrugBank XML data download."""
        #
        try:
            drugBankTargetMapPath = self.__getTargetDrugMapPath()
            self.__version = self.__dbP.getVersion()
            #
            if addTaxonomy:
                umP = UniProtIdMappingProvider(self.__cachePath)
                umP.reload(useCache=True)
            #
            oD = {}
            uD = {}
            taxonD = {}
            #
            dbIdL = self.__dbP.getDrugbankIds()
            logger.info("Building resource file for (%d) DrugBank entries", len(dbIdL))
            for dbId in dbIdL:
                tiDL = self.__dbP.getFeature(dbId, "target_interactions")
                for tiD in tiDL:
                    dD = {}
                    dD["drugbank_id"] = dbId
                    dD["name"] = self.__dbP.getFeature(dbId, "name")
                    dD["description"] = self.__dbP.getFeature(dbId, "description")
                    dD["moa"] = self.__dbP.getFeature(dbId, "mechanism-of-action")
                    dD["pharmacology"] = self.__dbP.getFeature(dbId, "pharmacodynamics")
                    dD["inchi_key"] = self.__dbP.getFeature(dbId, "inchikey")
                    dD["smiles"] = self.__dbP.getFeature(dbId, "smiles")
                    dD["pubmed_ids"] = tiD["articles"]
                    tS = tiD["amino-acid-sequence"]
                    if not tS:
                        logger.debug("Skipping entry target %r", tiD)
                        continue
                    sequence = "".join(tS.split("\n")[1:])
                    unpId = tiD["uniprot_ids"]
                    if "," in unpId or ";" in unpId or isinstance(unpId, list):
                        logger.warning("Bad uniprot id %r", unpId)
                    dD["target_name"] = tiD["name"]
                    cD = {"sequence": sequence, "uniprotId": unpId}
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
                        taxonD["%s\t%s" % (seqId, taxId)] = True
                    #
                    logger.debug("%r dD %r", dbId, dD)
                    uD.setdefault(unpId, []).append(dD)

            ok1 = self.__mU.doExport(fastaPath, oD, fmt="fasta", makeComment=True)
            tS = datetime.datetime.now().isoformat()
            vS = self.__version
            ok2 = self.__mU.doExport(drugBankTargetMapPath, {"version": vS, "created": tS, "cofactors": uD}, fmt="json", indent=3)
            ok3 = True
            if addTaxonomy:
                ok3 = self.__mU.doExport(taxonPath, list(taxonD.keys()), fmt="list")
            return ok1 & ok2 & ok3
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def __decodeComment(self, comment, separator="|"):
        dD = {}
        try:
            ti = iter(comment.split(separator))
            dD = {tup[1]: tup[0] for tup in zip(ti, ti)}
        except Exception:
            pass
        return dD
