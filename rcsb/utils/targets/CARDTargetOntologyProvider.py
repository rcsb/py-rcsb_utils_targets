##
#  File:           CARDTargetOntologyProvider.py
#  Date:           14-Mar-2023 dwp
#
#  Updates:
#   11-Apr-2023 dwp  Fix issue with lineage tree building--handle cases with two parents at same depth;
#                    Add treeNodeList building and exporting
#   27-Apr-2023 dwp  Update tree node list generation
#    2-May-2023 dwp  Remove depth field from lineage tree
##
"""
Accessors for CARD ontologies.

"""

import datetime
import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class CARDTargetOntologyProvider:
    """Accessors for CARD ontologies."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "CARD-ontology")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__oD, self.__tnL, self.__version = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=500):
        if self.__oD and len(self.__oD) > minCount:
            return True
        else:
            return False

    def getAssignmentVersion(self):
        return self.__version

    def getLineage(self, aroId):
        """Return the lineage (all parents + the requested aroId) of the given aroId.

        Args:
            aroId (str): ARO ID in the form of "ARO:3001059"

        Returns:
            list: list of dictionaries containing the "id" and "name" of each parent
        """
        return self.__oD.get(aroId, [])

    def getTreeNodeList(self):
        return self.__tnL

    def __reload(self, dirPath, **kwargs):
        oD = None
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
            oD = qD["data"]
        else:
            logger.info("Fetching url %s path %s", ontologyDumpUrl, ontologyDumpPath)
            ok = fU.get(ontologyDumpUrl, ontologyDumpPath)
            fU.mkdir(ontologyDumpDirPath)
            fU.uncompress(ontologyDumpPath, outputDir=ontologyDumpDirPath)
            fU.unbundleTarfile(os.path.join(ontologyDumpDirPath, ontologyDumpFileName[:-4]), dirPath=ontologyDumpDirPath)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            oD, tnL, version = self.__parseOntologyData(os.path.join(ontologyDumpDirPath, "aro.obo"))
            #
            tS = datetime.datetime.now().isoformat()
            qD = {"version": version, "created": tS, "data": oD}
            oD = qD["data"]
            ok = self.__mU.doExport(ontologyDataPath, qD, fmt="json", indent=3)
            logger.info("Export CARD ontology data (%d) status %r", len(oD), ok)

        return oD, tnL, version

    def __parseOntologyData(self, filePath):
        """Parse CARD ontology data

        Args:
            filePath (str): card ontology data file (aro.obo)

        Returns:
            (dict, string): dictionary of all ARO IDs and their ancestor lineages, ontology version string
        """
        try:
            cpD, tnL = None, None
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

            cpD, tnL = self.__buildLineageTree(parentChildList, idNameD)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return cpD, tnL, version

    def __buildLineageTree(self, parentChildTupleList, idNameMapD):
        """Build a lineage tree containing all children as keys and a list of all possible parents as the values.

        Args:
            parentChildTupleList (list): list of (parent, child) tuples (e.g., [('ARO:1000003', 'ARO:0000000'), ('ARO:1000003', 'ARO:0000001')])
            idNameMapD (dict): dictionary mapping of ARO IDs to their corresponding name

        Returns:
            dict: dictionary containing all children as keys and all possible parents as values
                  (including the child itself, but excluding the top-level parent 'ARO:1000001')
            list: list of all nodes (as dicts) in the tree with their immediate parents only (for building tree in browser)
        """
        # create a dictionary to store the parents of each child
        childToParentD = {}
        for parent, child in parentChildTupleList:
            if child not in childToParentD and child != "ARO:1000001":
                childToParentD[child] = []
            if parent != "ARO:1000001":  # Exclude the top-level "ARO:1000001"
                childToParentD[child].append(parent)

        treeNodeL = self.__exportTreeNodeList(childToParentD, idNameMapD)

        # create a dictionary to store the ancestors of each child
        lineageD = {}
        for child in childToParentD:
            lineageD[child] = [{"id": child, "name": idNameMapD[child]}]  # Add the child to its own ancestry list
            stack = [child]
            while stack:
                node = stack.pop()
                if node in childToParentD:
                    for parent in childToParentD[node]:
                        if parent not in [d["id"] for d in lineageD[child]]:
                            lineageD[child].append({"id": parent, "name": idNameMapD[parent]})
                            stack.append(parent)

        return lineageD, treeNodeL

    def __exportTreeNodeList(self, childToParentD, idNameMapD):
        """Create tree node list in the format of:

        {'id': 'ARO:1000003', 'name': 'antibiotic molecule'}
        {'id': 'ARO:0000041', 'name': 'bacitracin', 'parents': ['ARO:3000053', 'ARO:3000707']}
        {'id': 'ARO:0000039', 'name': 'spectinomycin', 'parents': ['ARO:0000016']}
        """
        #
        dL = []
        for child, parentL in childToParentD.items():
            if parentL:
                tD = {"id": child, "name": idNameMapD[child], "parents": parentL}
            else:
                tD = {"id": child, "name": idNameMapD[child]}
            dL.append(tD)

        return dL
