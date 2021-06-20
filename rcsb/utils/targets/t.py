class ChEMBLTargetActivityWorker(object):
    """A skeleton worker class that implements the interface expected by the multiprocessing module
    for fetching ChEMBL activity data --
    """

    def __init__(self, **kwargs):
        _ = kwargs

    def fetchActivity(self, dataList, procName, optionsD, workingDir):
        """Fetch ChEMBL activity for the input ChEMBL target identifier list"""
        _ = workingDir
        successList = []
        failList = []
        retList = []
        diagList = []
        #
        try:
            chunkSize = optionsD.get("chunkSize", 50)
            atL = optionsD.get("attributeList", [])

            for ii in range(0, len(dataList), chunkSize):
                logger.info("Begin chunk at ii %d/%d", ii, len(dataList))
                act = new_client.activity  # pylint: disable=no-member
                act.set_format("json")
                actDL = (
                    act.filter(target_chembl_id__in=dataList[ii : ii + chunkSize]).filter(standard_type__in=["IC50", "Ki", "EC50", "Kd"]).filter(standard_value__isnull=False).only(atL)
                )
                logger.info("Results (%d)", len(actDL))
                if actDL:
                    retList.extend(self.__activitySelect(atL, actD))
            successList = sorted(set(dataList) - set(failList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s built target interactions for %d/%d entries failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))
        #
        return successList, retList, diagList
