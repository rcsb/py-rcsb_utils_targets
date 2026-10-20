##
#  File:           TargetCofactorDbProvider.py
#  Date:           20-Aug-2024 dwp
#
#  Updated:
#  23-Aug-2024 dwp Use context manager for opening connections to ensure they close properly
#  26-Aug-2024 dwp Revert back to persistent client connections
#  14-Jul-2026 dwp Stop forcing upper case entity IDs when loading and fetching cofactor data, to handle short or
#                  extended PDB IDs
#  20-Oct-2026 dwp Open the persistent client lazily, once per process, so that no client is carried across a
#                  fork into multiprocessing workers; use a short-lived connection for document counts
#
##
"""
Provider/accessors for target cofactors data to/from MongoDB.
"""

import logging
import os
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
        self.__chunkSize = kwargs.get("chunkSize", 50)
        self.__resourceName = "MONGO_DB"
        #
        self.__databaseName = kwargs.get("databaseName", "cofactor_exdb")
        self.__collectionName = kwargs.get("collectionName", self.__cofactorResourceName)
        #
        # The persistent client is opened on first use in each process (see __getCollection())
        self.__client = None
        self.__clientPid = None

    def __getCollection(self):
        """Return the cofactor collection, opening a persistent client on first use in the current process.

        Loaders create this provider in the parent process and then fork worker processes (MultiProcUtil), and
        a MongoClient shouldn't be carried across a fork. Tying the client to the process ID gives each worker
        its own client, opened once and reused for all of that worker's lookups.
        """
        pid = os.getpid()
        if self.__client is None or self.__clientPid != pid:
            # Drop (don't close) a client inherited from a parent process, since closing it here would act on the
            # parent's copy of the client state
            self.__client = None
            conn = Connection(cfgOb=self.__cfgOb, resourceName=self.__resourceName)
            conn.openConnection()
            client = conn.getClientConnection()
            if client is None:
                raise ConnectionError("Unable to connect to resource %s" % self.__resourceName)
            self.__client = client
            self.__clientPid = pid
        return self.__client[self.__databaseName][self.__collectionName]

    def cofactorDbCount(self):
        """Count the number of documents in the cofactor collection.

        Uses a short-lived connection, so that checking the collection (e.g., in testCache()) leaves no open client
        behind to be inherited by forked worker processes.

        Returns:
            int: document count, or -1 if the collection can't be reached
        """
        count = -1
        try:
            with Connection(cfgOb=self.__cfgOb, resourceName=self.__resourceName) as client:
                if client is None:
                    raise ConnectionError("Unable to connect to resource %s" % self.__resourceName)
                count = client[self.__databaseName][self.__collectionName].count_documents({})
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return count

    def fetchCofactorData(self, rcsbEntityId, dataFieldName="rcsb_cofactors"):
        """Fetches the 'rcsb_cofactors' field for a document with the given rcsb_id.

        :param rcsb_id: The rcsb_id of the document to fetch.
        :return: The 'rcsb_cofactors' field if present, otherwise None.
        """
        data = []
        try:
            document = self.__getCollection().find_one({"rcsb_id": rcsbEntityId})
            if document and dataFieldName in document:
                data = document[dataFieldName]
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return data

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
        dL = [{"rcsb_id": k, "rcsb_cofactors": vL} for k, vL in cofactorDataDict.items()]
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
