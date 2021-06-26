##
#  File:           ChEMBLTargetProvider.py
#  Date:           9-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for ChEMBL target assignments.

"""

import logging
import os.path
import time

from chembl_webresource_client.new_client import new_client

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seq.UniProtIdMappingProvider import UniProtIdMappingProvider


logger = logging.getLogger(__name__)


class ChEMBLTargetProvider:
    """Accessors for ChEMBL target assignments."""

    def __init__(self, cachePath, useCache, **kwargs):
        #
        self.__cachePath = cachePath
        self.__dirPath = os.path.join(self.__cachePath, "ChEMBL-targets")
        baseVersion = 28
        self.__version = baseVersion
        self.__mapD = self.__reload(self.__dirPath, baseVersion, useCache, **kwargs)
        #

    def testCache(self, minCount=0):
        return self.__mapD and len(self.__mapD) > minCount

    def getAssignmentVersion(self):
        return self.__version

    def getTargetDescription(self, unpId):
        try:
            return self.__mapD[unpId][1]
        except Exception:
            logger.error("Missing description for %r", unpId)
            return None

    def getTargetChEMBLId(self, unpId):
        try:
            return self.__mapD[unpId][0]
        except Exception:
            return None

    def getTargetDataPath(self):
        return os.path.join(self.__dirPath, "chembl-target-data.json")

    def __reload(self, dirPath, baseVersion, useCache, **kwargs):
        startTime = time.time()
        mU = MarshalUtil(workPath=dirPath)
        chemblDbUrl = kwargs.get("ChEMBLDbUrl", "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/")
        ok = False
        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        # ChEMBL current version <baseVersion>,...
        # template:  chembl_<baseVersion>.fa.gz
        #
        targetFileName = "chembl_" + str(baseVersion) + ".fa.gz"
        mappingFileName = "chembl_uniprot_mapping.txt"
        #
        chemblTargetPath = os.path.join(dirPath, targetFileName)
        chemblMappingPath = os.path.join(dirPath, mappingFileName)
        mappingFilePath = os.path.join(dirPath, "chembl_uniprot_mapping.json")
        #
        mapD = {}
        if useCache and fU.exists(mappingFilePath):
            logger.info("useCache %r using %r and %r and %r", useCache, chemblTargetPath, chemblMappingPath, mappingFilePath)
            mapD = mU.doImport(mappingFilePath, fmt="json")
        else:
            # Get the ChEMBL UniProt mapping file
            url = os.path.join(chemblDbUrl, mappingFileName)
            ok = fU.get(url, chemblMappingPath)
            logger.info("Fetched %r url %s path %s", ok, url, chemblMappingPath)
            logger.info("Reading ChEMBL mapping file path %s", mappingFilePath)
            rowL = mU.doImport(chemblMappingPath, fmt="tdd", rowFormat="list")
            for row in rowL:
                mapD[row[0]] = (row[1], row[2], row[3])
            ok = mU.doExport(mappingFilePath, mapD, fmt="json")
            logger.info("Processed mapping path %s (%d) %r", mappingFilePath, len(mapD), ok)
            #
            # Get the target FASTA files --
            for vers in range(baseVersion, baseVersion + 10):
                logger.info("Now fetching version %r", vers)
                self.__version = vers
                targetFileName = "chembl_" + str(vers) + ".fa.gz"
                chemblTargetPath = os.path.join(dirPath, "chembl_targets_raw.fa.gz")
                url = os.path.join(chemblDbUrl, targetFileName)
                ok = fU.get(url, chemblTargetPath)
                logger.info("Fetched %r url %s path %s", ok, url, chemblTargetPath)
                if ok:
                    break
        #
        logger.info("Completed reload at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        #
        return mapD

    def exportFasta(self, fastaPath, taxonPath, addTaxonomy=False):
        ok = self.__parseFasta(fastaPath, taxonPath, self.__cachePath, self.__dirPath, addTaxonomy=addTaxonomy)
        return ok

    def __parseFasta(self, fastaPath, taxonPath, cachePath, dirPath, addTaxonomy=False):
        # input paths
        chemblTargetRawPath = os.path.join(dirPath, "chembl_targets_raw.fa.gz")
        mU = MarshalUtil(workPath=cachePath)
        oD = {}
        uD = {}
        missTax = 0
        taxonL = []
        try:
            if addTaxonomy:
                umP = UniProtIdMappingProvider(cachePath)
                umP.reload(useCache=True)
            #
            fD = mU.doImport(chemblTargetRawPath, fmt="fasta", commentStyle="default")
            #
            for seqId, sD in fD.items():
                chemblId = seqId.strip().split(" ")[0].strip()
                unpId = seqId[seqId.find("[") + 1 : seqId.find("]")]
                seq = sD["sequence"]
                cD = {"sequence": seq, "uniprotId": unpId, "chemblId": chemblId}
                if addTaxonomy:
                    taxId = umP.getMappedId(unpId, mapName="NCBI-taxon")
                    cD["taxId"] = taxId if taxId else -1
                    if not taxId:
                        missTax += 1
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
                #
                uD.setdefault(unpId, []).append(chemblId)
            #
            ok1 = mU.doExport(fastaPath, oD, fmt="fasta", makeComment=True)
            ok3 = True
            if addTaxonomy:
                ok3 = mU.doExport(taxonPath, taxonL, fmt="list")
            return ok1 & ok3
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return False

    def getActivityData(self, targetChEMBLIdList):
        """Get cofactor activity data for the input ChEMBL target list.

        Args:
            targetChEMBLIdList (list): list of ChEMBL target identifiers

        Returns:
          (dict, dict):  {targetChEMBId: {activity data}}, {moleculeChEMBId: {activity data}}

        """
        targetD = {}
        chunkSize = 50
        numTargets = len(targetChEMBLIdList)
        try:
            for ii in range(0, len(targetChEMBLIdList), chunkSize):
                logger.info("Begin chunk ii %d/%d", ii, numTargets)
                act = new_client.activity  # pylint: disable=no-member
                act.set_format("json")
                actDL = (
                    act.filter(target_chembl_id__in=targetChEMBLIdList[ii : ii + chunkSize]).filter(standard_type__in=["IC50", "Ki", "EC50", "Kd"]).filter(standard_value__isnull=False)
                )
                if actDL:
                    for actD in actDL:
                        targetD.setdefault(actD["target_chembl_id"], []).append(self.__activitySelect(actD))
                logger.info("End chunk completed (%d)", ii)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return targetD

    def __activitySelect(self, aD):
        atL = [
            "assay_chembl_id",
            "assay_description",
            "assay_type",
            "canonical_smiles",
            "ligand_efficiency",
            "molecule_chembl_id",
            "parent_molecule_chembl_id",
            "pchembl_value",
            "standard_relation",
            "standard_type",
            "standard_units",
            "standard_value",
            "target_chembl_id",
        ]
        return {at: aD[at] if at in aD else None for at in atL}

    def getMechanismData(self, targetChEMBLIdList):
        """Get mechanism data for the input ChEMBL target list.

        Args:
            targetChEMBLIdList (list): list of ChEMBL target identifiers

        Returns:
          (dict):  dictionary  {ChEMBId: {mechanism data}}

        """
        oD = {}
        chunkSize = 50
        try:
            for ii in range(0, len(targetChEMBLIdList), chunkSize):
                mch = new_client.mechanism  # pylint: disable=no-member
                mch.set_format("json")
                mDL = mch.filter(target_chembl_id__in=targetChEMBLIdList[ii : ii + chunkSize])
                if mDL:
                    logger.info("mDL (%d)", len(mDL))
                    for mD in mDL:
                        oD.setdefault(mD["target_chembl_id"], []).append(mD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD
