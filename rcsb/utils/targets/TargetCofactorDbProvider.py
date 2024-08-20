##
#  File:           TargetCofactorDbProvider.py
#  Date:           20-Aug-2024 dwp
#
#  Updated:
#
##
"""
Provider/accessors for target cofactors data to/from MongoDB.
"""

import logging
import time

from rcsb.db.mongo.Connection import Connection
from rcsb.db.mongo.DocumentLoader import DocumentLoader

logger = logging.getLogger(__name__)


class TargetCofactorDbProvider:
    def __init__(self, cachePath, cfgOb=None, cofactorResourceName=None, **kwargs):
        """
        Provider class for loading and fetching cofactor data to/from MongoDB.
        """
        self.__cachePath = cachePath
        self.__cfgOb = cfgOb
        self.__cofactorResourceName = cofactorResourceName  # "chembl", "pharos", or "drugbank"
        self.__numProc = kwargs.get("numProc", 6)
        self.__chunkSize = kwargs.get("chunkSize", 10)
        self.__resourceName = "MONGO_DB"
        #
        self.__databaseName = kwargs.get("databaseName", "cofactor_exdb")
        self.__collectionName = kwargs.get("collectionName", self.__cofactorResourceName)
        #
        conn = Connection(cfgOb=self.__cfgOb, resourceName=self.__resourceName)
        conn.openConnection()
        self.__client = conn.getClientConnection()
        self.__db = self.__client[self.__databaseName]
        self.__collection = self.__db[self.__collectionName]

    def cofactorDbCount(self):
        """Count the number of documents in the cofactor collection.
        """
        count = 0
        try:
            count = self.__collection.count_documents({})
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return count

    def fetchCofactorData(self, rcsbEntityId, dataFieldName="rcsb_cofactors"):
        """Fetches the 'rcsb_cofactors' field for a document with the given rcsb_id.

        :param rcsb_id: The rcsb_id of the document to fetch.
        :return: The 'rcsb_cofactors' field if present, otherwise None.
        """
        document = self.__collection.find_one({"rcsb_id": rcsbEntityId.upper()})
        if document and dataFieldName in document:
            return document[dataFieldName]
        else:
            return []

    def loadCofactorData(self, cofactorResourceName, cofactorProvider, **kwargs):
        """Load cofactor data for the input data resource.

        Args:
            cofactorResourceName (str): target cofactor resource name (e.g., "pharos", "chembl", or "drugbank).
            cofactorProvider (obj): target cofactor resource provider instance.

        Returns:
            bool: True for success or False otherwise
        """
        if cofactorResourceName not in ["chembl", "pharos", "drugbank"]:
            logger.error("Unsupported cofactor resource %r", cofactorResourceName)
            return False
        #
        loadType = kwargs.get("loadType", "full")
        maxStepLength = kwargs.get("maxStepLength", 10000)
        documentLimit = kwargs.get("documentLimit", None)
        verbose = kwargs.get("verbose", False)
        readBackCheck = kwargs.get("readBackCheck", False)
        #
        dl = DocumentLoader(
            self.__cfgOb,
            self.__cachePath,
            self.__resourceName,
            numProc=self.__numProc,
            chunkSize=self.__chunkSize,
            maxStepLength=maxStepLength,
            documentLimit=documentLimit,
            verbose=verbose,
            readBackCheck=readBackCheck,
        )
        #
        cofactorDataDict = cofactorProvider.getCofactorDataDict()
        dL = []
        dL = [{"rcsb_id": k.upper(), "rcsb_cofactors": vL} for k, vL in cofactorDataDict.items()]
        #
        startTime = time.time()
        ok = dl.load(self.__databaseName, self.__collectionName, loadType=loadType, documentList=dL, indexAttributeList=["rcsb_id"], keyNames=None, schemaLevel=None)
        logger.info(
            "Completed loading cofactor data for %s (%r) to databaseName %s collectionName %s (status %r) at %s (%.4f seconds)",
            cofactorResourceName,
            len(dL),
            self.__databaseName,
            self.__collectionName,
            ok,
            time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
            time.time() - startTime,
        )
        #
        return ok
