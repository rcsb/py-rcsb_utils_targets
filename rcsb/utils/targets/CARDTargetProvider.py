##
#  File:           CARDTargetProvider.py
#  Date:           27-Nov-2020 jdw
#
#  Updated:
#    6-Mar-2023 dwp  Add filterForHomologs() method to filter CARD data for only protein homolog models;
#                    Add ontology data and parsing methods for getting lineage
##
"""
Accessors for CARD target assignments.

"""

import datetime
import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class CARDTargetProvider:
    """Accessors for CARD target assignments."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "CARD-targets")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__oD, self.__version = self.__reload(self.__dirPath, **kwargs)
        self.__cpD, self.__ontologyVersion = self.__reloadOntology(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=3000):
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def hasModel(self, modelId):
        return modelId in self.__oD

    def getModelKey(self, modelId, key):
        try:
            return self.__oD[modelId][key]
        except Exception:
            return None

    def getAssignmentVersion(self):
        return self.__version

    def getOntologyVersion(self):
        return self.__ontologyVersion

    def getLineage(self, aroId):
        """Return the lineage (all parents + the requested aroId) of the given aroId.

        Args:
            aroId (str): ARO ID in the form of "ARO:3001059"

        Returns:
            list: list of dictionaries containing the "id" and "name" of each parent
        """
        return self.__cpD.get(aroId, [])

    def getTargetDataPath(self):
        return os.path.join(self.__dirPath, "card-target-data.json")

    def getCofactorDataPath(self):
        return None

    def __reload(self, dirPath, **kwargs):
        oD = None
        version = None
        startTime = time.time()
        useCache = kwargs.get("useCache", True)
        #
        ok = False
        fU = FileUtil()
        #
        # CARDDumpUrl = kwargs.get("CARDDumpUrl", "https://card.mcmaster.ca/latest/data/broadstreet-v3.1.0.tar.bz2")
        cardDumpUrl = kwargs.get("CARDDumpUrl", "https://card.mcmaster.ca/latest/data")
        cardDumpFileName = "card-data.tar.bz2"
        cardDumpPath = os.path.join(dirPath, cardDumpFileName)
        cardDumpDirPath = os.path.join(dirPath, "dump-targets")
        #
        fU.mkdir(dirPath)
        cardDataPath = os.path.join(dirPath, "card-select-data.json")
        #
        # ---
        # Load CARD data
        logger.info("useCache %r CARDDumpPath %r", useCache, cardDumpPath)
        if useCache and self.__mU.exists(cardDataPath):
            qD = self.__mU.doImport(cardDataPath, fmt="json")
            version = qD["version"]
            oD = qD["data"]
        else:
            logger.info("Fetching url %s path %s", cardDumpUrl, cardDumpPath)
            ok = fU.get(cardDumpUrl, cardDumpPath)
            fU.mkdir(cardDumpDirPath)
            fU.uncompress(cardDumpPath, outputDir=cardDumpDirPath)
            fU.unbundleTarfile(os.path.join(cardDumpDirPath, cardDumpFileName[:-4]), dirPath=cardDumpDirPath)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            oD, version = self.__parseCardData(os.path.join(cardDumpDirPath, "card.json"))
            #
            # Filter for only protein homolog models (this can be removed if/when we decide to include all types of models)
            oD = self.__filterForHomologs(os.path.join(cardDumpDirPath, "protein_fasta_protein_homolog_model.fasta"), oD)
            #
            tS = datetime.datetime.now().isoformat()
            qD = {"version": version, "created": tS, "data": oD}
            oD = qD["data"]
            ok = self.__mU.doExport(cardDataPath, qD, fmt="json", indent=3)
            logger.info("Export CARD data (%d) status %r", len(oD), ok)
        # ---
        return oD, version

    def __reloadOntology(self, dirPath, **kwargs):
        cpD = None
        version = None
        startTime = time.time()
        useCache = kwargs.get("useCache", True)
        #
        ok = False
        fU = FileUtil()
        #
        ontologyDumpUrl = kwargs.get("OntologyDumpUrl", "https://card.mcmaster.ca/latest/ontology")
        ontologyDumpFileName = "card-ontology.tar.bz2"
        ontologyDumpPath = os.path.join(dirPath, ontologyDumpFileName)
        ontologyDumpDirPath = os.path.join(dirPath, "dump-ontology")
        #
        fU.mkdir(dirPath)
        ontologyDataPath = os.path.join(dirPath, "card-ontology-data.json")
        #
        # ---
        # Load Ontology data
        logger.info("useCache %r ontologyDumpPath %r", useCache, ontologyDumpPath)
        if useCache and self.__mU.exists(ontologyDataPath):
            qD = self.__mU.doImport(ontologyDataPath, fmt="json")
            version = qD["version"]
            cpD = qD["data"]
        else:
            logger.info("Fetching url %s path %s", ontologyDumpUrl, ontologyDumpPath)
            ok = fU.get(ontologyDumpUrl, ontologyDumpPath)
            fU.mkdir(ontologyDumpDirPath)
            fU.uncompress(ontologyDumpPath, outputDir=ontologyDumpDirPath)
            fU.unbundleTarfile(os.path.join(ontologyDumpDirPath, ontologyDumpFileName[:-4]), dirPath=ontologyDumpDirPath)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            cpD, version = self.__parseOntologyData(os.path.join(ontologyDumpDirPath, "aro.obo"))
            #
            tS = datetime.datetime.now().isoformat()
            qD = {"version": version, "created": tS, "data": cpD}
            cpD = qD["data"]
            ok = self.__mU.doExport(ontologyDataPath, qD, fmt="json", indent=3)
            logger.info("Export CARD ontology data (%d) status %r", len(cpD), ok)

        return cpD, version

    def exportCardFasta(self, fastaPath, taxonPath):
        ok = self.__exportCardFasta(fastaPath, taxonPath, self.__oD)
        return ok

    def __exportCardFasta(self, fastaPath, taxonPath, cardD):
        """Export a CARD sequence target fasta file

        Args:
            fastaPath (str): fasta output file path
            cardD (dict): card selected data dictionary

        Returns:
            (bool): True for success or False otherwise
        """
        sD = {}
        taxonL = []
        try:
            for modelId, tD in cardD.items():
                modelBitScore = None
                # aroAcc = tD["accession"]
                aroId = tD["cvTermId"]
                if "sequences" not in tD:
                    continue
                modelBitScore = tD["modelBitScore"] if "modelBitScore" in tD else None
                for qD in tD["sequences"]:
                    sId = qD["seqId"]
                    seq = qD["sequence"]
                    taxId = qD["taxId"]
                    cD = {"sequence": seq, "modelId": modelId, "aroId": aroId, "seqId": sId, "taxId": taxId}
                    cD["bitScore"] = modelBitScore if modelBitScore else "-1.0"
                    #
                    cId = ""
                    cL = []
                    for k, v in cD.items():
                        if k in ["sequence"]:
                            continue
                        cL.append(str(v))
                        cL.append(str(k))
                    cId = "|".join(cL)
                    sD[cId] = cD
                    taxonL.append("%s\t%s" % (cId, taxId))

            ok = self.__mU.doExport(fastaPath, sD, fmt="fasta", makeComment=True)
            logger.info("Export CARD fasta (%d) status %r", len(sD), ok)
            ok = self.__mU.doExport(taxonPath, taxonL, fmt="list")
            logger.info("Export Taxon (%d) status %r", len(taxonL), ok)
        except Exception as e:
            logger.exception("Failing for model %r tD %r with %s", modelId, tD, str(e))
        return ok

    def __parseCardData(self, filePath):
        """Parse CARD target data

        Args:
            filePath (str): card json data file

        Returns:
            (dict, string): card selected data dictionary, card version string
        """
        try:
            oD = {}
            version = None
            cD = self.__mU.doImport(filePath, fmt="json")
            logger.info("CARD model count (%d)", len(cD))
            for modelId, mD in cD.items():
                if modelId.startswith("_"):
                    if modelId == "_version":
                        version = mD
                    continue
                oD[modelId] = {}  # modelId = '1028'
                for kTup in [
                    ("ARO_accession", "accession"),  # 'ARO_accession', '3001059'
                    ("ARO_id", "cvTermId"),  # 'ARO_id', '37439'
                    ("ARO_name", "name"),  # 'ARO_name', 'SHV-1'
                    ("ARO_description", "descr"),  # 'SHV-1 is a broad-spectrum beta-lactamase...
                    ("model_name", "modelName"),  # 'model_name', 'SHV-1'
                    ("model_type", "modelType"),  # 'model_type', 'protein homolog model'
                ]:
                    if kTup[0] in mD:
                        oD[modelId][kTup[1]] = mD[kTup[0]]
                # Add in category details (Family, Drug Class(es), and Resistance mechanism)
                acD = {}
                if "ARO_category" in mD:
                    aroCategoryD = mD["ARO_category"]  # 'ARO_category', large dict containing each type of 'category_aro_class_name' ("AMR Gene Family", "Drug Class", ...)
                    for cvId, catD in aroCategoryD.items():
                        if all([k in catD for k in ["category_aro_class_name", "category_aro_accession", "category_aro_name", "category_aro_description"]]):
                            if catD["category_aro_class_name"] == "AMR Gene Family":
                                acD["familyCvTermId"] = cvId
                                acD["familyAccession"] = catD["category_aro_accession"]
                                acD["familyName"] = catD["category_aro_name"]
                                acD["familyDescription"] = catD["category_aro_description"]
                            if catD["category_aro_class_name"] == "Drug Class":
                                acD.setdefault("drugClasses", []).append(catD["category_aro_name"])
                            if catD["category_aro_class_name"] == "Resistance Mechanism":
                                acD["resistanceMechanism"] = catD["category_aro_name"]
                oD[modelId].update(acD)
                #
                try:
                    if "model_sequences" in mD:
                        for seqId, tD in mD["model_sequences"]["sequence"].items():
                            oD[modelId].setdefault("sequences", []).append(
                                {"seqId": seqId, "sequence": tD["protein_sequence"]["sequence"], "taxId": tD["NCBI_taxonomy"]["NCBI_taxonomy_id"]}
                            )
                except Exception as e:
                    logger.exception("Failing with %s", str(e))

                try:
                    if "model_param" in mD and "blastp_bit_score" in mD["model_param"] and "param_value" in mD["model_param"]["blastp_bit_score"]:
                        oD[modelId]["modelBitScore"] = mD["model_param"]["blastp_bit_score"]["param_value"]

                except Exception as e:
                    logger.exception("Failing with %s", str(e))

        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD, version

    def __filterForHomologs(self, filePath, targetD):
        """Filter the CARD target data to select only for protein homologs

        Args:
            filePath (str): path to CARD fasta protein homoolog model file (protein_fasta_protein_homolog_model.fasta; from source)
            targetD (dict): dictionary generated from __parseCardData (i.e., oD)

        Returns:
            (dict): filtered card selected data dictionary
        """

        with open(filePath, "r", encoding="utf-8") as f:
            data = f.readlines()

        aroL = [i.split("|")[2].strip("ARO:") for i in data if i.startswith(">")]

        filteredD = {}
        for k, v in targetD.items():
            if v["accession"] in aroL:
                filteredD.update({k: v})

        return filteredD

    def __parseOntologyData(self, filePath):
        """Parse CARD ontology data

        Args:
            filePath (str): card ontology data file (aro.obo)

        Returns:
            (dict, string): dictionary of all ARO IDs and their ancestor lineages, ontology version string
        """
        try:
            cpD = None
            version = None
            parentChildList = []
            idNameD = {}

            with open(filePath, "r", encoding="utf-8") as f:
                data = f.read()

            dL = data.split("[Term]")

            version = dL.pop(0).split("\n")[0].split("format-version: ")[1]

            for tS in dL:
                tsL = tS.strip().split("\n")
                childId, parentId, childName = None, None, None
                for item in tsL:
                    if item.startswith("id: "):
                        childId = item.split('id: ')[1]
                    if item.startswith("is_a: "):
                        parentId = item.split('is_a: ')[1].split(" !")[0]
                    if item.startswith("name: "):
                        childName = item.split('name: ')[1]
                    if parentId and childId and (parentId, childId) not in parentChildList:
                        parentChildList.append((parentId, childId))
                    if childName and childId not in idNameD:
                        idNameD[childId] = childName
                    if item.startswith("[Typedef]"):
                        break

            logger.info("Ontology parent-child pair count (%d)", len(parentChildList))

            cpD = self.__buildLineageTree(parentChildList, idNameD)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return cpD, version

    def __buildLineageTree(self, parentChildTupleList, idNameMapD):
        """Build a lineage tree containing all children as keys and a
        list of all possible parents as the values.

        Args:
            parentChildTupleList (list): list of (parent, child) tuples (e.g., [('ARO:1000003', 'ARO:0000000'), ('ARO:1000003', 'ARO:0000001')])
            idNameMapD (dict): dictionary mapping of ARO IDs to their corresponding name

        Returns:
            dict: dictionary containing all children as keys and all possible parents as values
                  (including the child itself, but excluding the top-level parent 'ARO:1000001')
        """
        # create a dictionary to store the parents of each child
        parents = {}
        for parent, child in parentChildTupleList:
            if child not in parents:
                parents[child] = []
            parents[child].append(parent)

        # create a dictionary to store the ancestors of each child
        lineageD = {}
        for child in parents.keys():
            lineageD[child] = []
            stack = [child]
            while stack:
                node = stack.pop()
                if node in parents:
                    for parent in parents[node]:
                        lineageD[child].append(parent)
                        stack.append(parent)

        # Go through the lineage dict and add the child to its own list and remove the top-level "ARO:1000001"
        # Also, update the list to contain both the parent ARO IDs and their associated names.
        for child, pL in lineageD.items():
            if child not in pL:
                npL = [child] + [p for p in pL if p != "ARO:1000001"]
                npL = npL.reverse()  # List oldest/broadest ancestor first (depth=1), youngest/most-specific child last (depth=n)
                npL = [{"id": id, "name": idNameMapD[id], "depth": ii + 1} for ii, id in enumerate(npL)]
                lineageD[child] = npL

        return lineageD
