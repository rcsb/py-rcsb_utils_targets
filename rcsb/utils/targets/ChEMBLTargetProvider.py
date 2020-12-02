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

    def __init__(self, **kwargs):
        #
        self.__reload(**kwargs)
        self.__parseFasta(**kwargs)
        #

    def testCache(self):
        return True

    def __reload(self, **kwargs):
        startTime = time.time()
        cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        chemblDbUrl = kwargs.get("ChEMBLDbUrl", "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/")
        dirPath = os.path.join(cachePath, "ChEMBL-targets")
        ok = False
        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        baseVersion = 27
        # ChEMBL current version 27,...
        # template:  chembl_27.fa.gz
        #
        targetFileName = "chembl_" + str(baseVersion) + ".fa.gz"
        mappingFileName = "chembl_uniprot_mapping.txt"
        #
        chemblTargetPath = os.path.join(dirPath, targetFileName)
        chemblMappingPath = os.path.join(dirPath, mappingFileName)
        #
        if useCache and fU.exists(chemblMappingPath):
            logger.info("useCache %r using %r and %r", useCache, chemblTargetPath, chemblMappingPath)
            ok = True
        else:
            url = os.path.join(chemblDbUrl, mappingFileName)
            ok = fU.get(url, chemblMappingPath)
            logger.info("Fetching url %s path %s", url, chemblMappingPath)
            #
            for vers in range(baseVersion, baseVersion + 10):
                targetFileName = "chembl_" + str(vers) + ".fa.gz"
                chemblTargetPath = os.path.join(dirPath, "chembl_targets_raw.fa.gz")
                url = os.path.join(chemblDbUrl, targetFileName)
                ok = fU.get(url, chemblTargetPath)
                logger.info("Fetching url %s path %s", url, chemblTargetPath)
                if ok:
                    break
            #
            logger.info("Completed fetches at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            #
        return ok

    def __parseFasta(self, **kwargs):
        cachePath = kwargs.get("cachePath", ".")
        dirPath = os.path.join(cachePath, "ChEMBL-targets")
        #
        chemblTargetRawPath = os.path.join(dirPath, "chembl_targets_raw.fa.gz")
        mappingFilePath = os.path.join(dirPath, "chembl_uniprot_mapping.txt")
        #
        chemblTargetPath = os.path.join(dirPath, "chembl_targets.fa")
        chemblTargetUniprotMapPath = os.path.join(dirPath, "chembl_uniprot_target_mapping.json")
        #
        chemblTargetUniprotRawMap = os.path.join(dirPath, "chembl_uniprot_target_mapping_raw.json")
        #
        addTaxonomy = kwargs.get("addTaxonomy", False)
        mU = MarshalUtil(workPath=cachePath)
        #
        mapD = {}
        rowL = mU.doImport(mappingFilePath, fmt="tdd", rowFormat="list")
        for row in rowL:
            mapD.setdefault(row[0], []).append((row[1], row[3]))
        ok = mU.doExport(chemblTargetUniprotRawMap, mapD, fmt="json")
        logger.info("ChEMBL Raw mapping path %s (%d) %r", chemblTargetUniprotRawMap, len(mapD), ok)
        #
        oD = {}
        uD = {}
        missTax = 0
        try:
            if addTaxonomy:
                umP = UniProtIdMappingProvider(cachePath=cachePath, useCache=True)
            #
            fD = mU.doImport(chemblTargetRawPath, fmt="fasta", commentStyle="default")
            #
            for seqId, sD in fD.items():
                chemblId = seqId.strip().split(" ")[0].strip()
                unpId = seqId[seqId.find("[") + 1 : seqId.find("]")]
                seq = sD["sequence"]
                oD[chemblId] = {"sequence": seq, "uniprotId": unpId, "chemblId": chemblId}
                if addTaxonomy:
                    taxId = umP.getMappedId(unpId, mapName="NCBI-taxon")
                    oD[chemblId]["taxId"] = taxId if taxId else -1
                    if not taxId:
                        missTax += 1
                #
                uD.setdefault(unpId, []).append(chemblId)
            #
            ok1 = mU.doExport(chemblTargetPath, oD, fmt="fasta", makeComment=True)
            ok2 = mU.doExport(chemblTargetUniprotMapPath, uD, fmt="json")
            logger.info("ChEMBL mapping path %s (%d) missing taxonomy count (%d)", chemblTargetUniprotMapPath, len(uD), missTax)
            return ok1 & ok2
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
        molD = {}
        chunkSize = 50
        try:
            for ii in range(0, len(targetChEMBLIdList), chunkSize):
                # for targetChEMBLId in targetChEMBLIdList:
                act = new_client.activity  # pylint: disable=no-member
                act.set_format("json")
                actDL = (
                    act.filter(target_chembl_id__in=targetChEMBLIdList[ii : ii + chunkSize])
                    .filter(standard_type__in=["IC50", "Ki", "EC50", "Kd"])
                    .filter(standard_value__isnull=False)
                )
                if actDL:
                    logger.info("actDL (%d)", len(actDL))
                    for actD in actDL:
                        targetD.setdefault(actD["target_chembl_id"], []).append(actD)
                        if "molecule_chembl_id" in actD:
                            molD.setdefault(actD["molecule_chembl_id"], []).append(actD)
                        if "parent_molecule_chembl_id" in actD and actD["molecule_chembl_id"] != actD["parent_molecule_chembl_id"]:
                            molD.setdefault(actD["parent_molecule_chembl_id"], []).append(actD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return targetD, molD

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
