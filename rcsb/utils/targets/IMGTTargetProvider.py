##
#  File:           IMGTTargetProvider.py
#  Date:           5-Jul-2021 jdw
#
#  Updated:
#
##
"""
Accessors for IMGT target annotations.
"""

import datetime
import glob
import gzip
import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.taxonomy.TaxonomyProvider import TaxonomyProvider

logger = logging.getLogger(__name__)


class IMGTTargetProvider(StashableBase):
    """Accessors for IMGT target annotations."""

    def __init__(self, cachePath, useCache, **kwargs):
        #
        self.__cachePath = cachePath
        self.__dirName = "IMGT-targets"
        imgtDumpUrl = kwargs.get("IMGTDumpUrl", None)
        super(IMGTTargetProvider, self).__init__(self.__cachePath, [self.__dirName])
        self.__dirPath = os.path.join(self.__cachePath, self.__dirName)
        self.__version = None
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__imgtD = self.__reload(self.__dirPath, useCache=useCache, imgtDumpUrl=imgtDumpUrl)
        #

    def testCache(self, minCount=1000):
        if self.__imgtD and "chains" in self.__imgtD and len(self.__imgtD["chains"]) > minCount:
            return True
        else:
            return False

    def getVersion(self):
        return self.__version

    def getChains(self):
        return self.__imgtD["chains"]

    def __reload(self, dirPath, useCache=False, imgtDumpUrl=None, testList=None, maxCount=None):
        imgtD = {}
        startTime = time.time()

        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        imgtDataPath = os.path.join(self.__dirPath, "imgt-data.json")
        #
        logger.info("useCache %r imgtFeaturePath %r", useCache, imgtDataPath)
        if useCache and self.__mU.exists(imgtDataPath):
            imgtD = self.__mU.doImport(imgtDataPath, fmt="json")
            self.__version = imgtD["version"]
        else:
            imgtDumpUrl = imgtDumpUrl if imgtDumpUrl else "http://www.imgt.org/download/3Dstructure-DB/IMGT3DFlatFiles.tgz"
            imgtReadmeUrl = "http://www.imgt.org/download/3Dstructure-DB/RELEASE"
            imgtDumpFileName = fU.getFileName(imgtDumpUrl)
            imgtDumpPath = os.path.join(dirPath, imgtDumpFileName)
            imgtReleasePath = os.path.join(dirPath, "IMGT-release.txt")
            _, fn = os.path.split(imgtDumpUrl)
            imgtFlatFilePath = os.path.join(self.__dirPath, fn[:-4])
            #
            logger.info("Fetching url %s path %s", imgtDumpUrl, imgtDumpPath)
            ok1 = fU.get(imgtDumpUrl, imgtDumpPath)
            ok2 = fU.get(imgtReadmeUrl, imgtReleasePath)
            fU.unbundleTarfile(imgtDumpPath, dirPath=dirPath)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok1 and ok2, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            # ---
            readmeLines = self.__mU.doImport(imgtReleasePath, fmt="list")
            self.__version = readmeLines[0].strip() if readmeLines else None
            logger.info("IMGT version %r", self.__version)
            # ---
            chainD, rawD = self.__imgtFlatFileProcessor(imgtFlatFilePath, maxCount=maxCount, testList=testList)
            # ---
            tS = datetime.datetime.now().isoformat()
            # vS = datetime.datetime.now().strftime("%Y-%m-%d")
            if testList:
                imgtD = {"version": self.__version, "date": tS, "chains": chainD, "raw": rawD}
            else:
                imgtD = {"version": self.__version, "date": tS, "chains": chainD}
            ok = self.__mU.doExport(imgtDataPath, imgtD, fmt="json", indent=3)
            logger.info("Completed flatfile prep (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return imgtD

    def exportFasta(self, withGaps=False):
        """
        Example:
            The IMGT/GENE-DB FASTA header contains 15 fields separated by '|':

            1. IMGT/LIGM-DB accession number(s)
            2. IMGT gene and allele name
            3. species (may be followed by an "_" and the name of the strain, breed or isolate, if defined)
            4. IMGT gene and allele functionality
            5. exon(s), region name(s), or extracted label(s)
            6. start and end positions in the IMGT/LIGM-DB accession number(s)
            7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
            8. codon start, or 'NR' (not relevant) for non coding labels
            9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
            10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
            11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
            12. number of amino acids (AA): this field indicates that the sequence is in amino acids
            13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
            14. partial (if it is)
            15. reverse complementary (if it is)

        """
        # --
        fU = FileUtil()
        fU.mkdir(self.__dirPath)
        if withGaps:
            imgtTargetUrl = "http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP"
        else:
            imgtTargetUrl = "http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP"
        imgtTargetFileName = fU.getFileName(imgtTargetUrl)
        rawFastaPath = os.path.join(self.__dirPath, imgtTargetFileName)
        # --
        logger.debug("Fetching url %s path %s", imgtTargetUrl, rawFastaPath)
        ok = fU.get(imgtTargetUrl, rawFastaPath)
        logger.info("Fetch status (%r) url %s path %s", ok, imgtTargetUrl, rawFastaPath)
        # --
        fastaPath = os.path.join(self.__dirPath, "imgt-reference.fa")
        taxonPath = os.path.join(self.__dirPath, "imgt-reference-taxon.tdd")
        tP = TaxonomyProvider(taxDirPath=self.__cachePath)

        rawQD = self.__mU.doImport(rawFastaPath, fmt="fasta", commentStyle="default")
        oD = {}
        taxonL = []
        for queryId, sD in rawQD.items():
            qL = queryId.split("|")
            tL = qL[2].split("_")
            taxName = tL[0]
            taxVar = tL[1].replace(" ", "_") if len(tL) > 1 else None
            taxId = tP.getTaxId(taxName)
            if taxId:
                tD = {"seqId": qL[0], "imgtGene": qL[1], "functionality": qL[3], "labels": qL[4], "taxId": taxId}
                if taxVar:
                    tD["taxVar"] = taxVar
                sD.update(tD)
            else:
                logger.info("Unknown taxonomy %r", queryId)
            sD["sequence"].replace(".", "-")
            seqId = ""
            cL = []
            for k, v in sD.items():
                if k in ["sequence"]:
                    continue
                cL.append(str(v))
                cL.append(str(k))
            seqId = "|".join(cL)
            oD[seqId] = sD
            taxonL.append("%s\t%s" % (seqId, taxId))
        #
        ok1 = self.__mU.doExport(taxonPath, taxonL, fmt="list")
        ok2 = self.__mU.doExport(fastaPath, oD, fmt="fasta", makeComment=True)
        return ok1 and ok2

    def __imgtFlatFileProcessor(self, flatFilePath, maxCount=None, testList=None):
        chainD = {}
        rawD = {}
        failures = []
        idList = []
        ic = 0
        filePattern = os.path.join(flatFilePath, "*.pdb.gz")
        logger.info("Collecting flat files with pattern %r", filePattern)
        for fp in glob.glob(filePattern):
            ic += 1
            if maxCount and ic > maxCount:
                break
            logger.debug("Processing file %r", fp)
            _, fn = os.path.split(fp)
            pdbId = fn[5:9].lower()
            if testList and pdbId not in testList:
                continue
            idList.append(pdbId)
            cD = {}
            tmpD = {}
            with gzip.open(fp, "rb") as ifh:
                try:
                    cD, tmpD = self.__imgtRemarkParser(pdbId, ifh)
                except Exception as e:
                    failures.append(pdbId)
                    logger.exception("Failing for %r with %s", pdbId, str(e))
                    continue
            # --
            chainD.update(cD)
            rawD[pdbId] = tmpD
        #
        logger.info("ID List (%d)", len(set(idList)))
        sL = list(rawD.keys())
        logger.info("Successes (%d) chains (%d)", len(sL), len(chainD))
        logger.info("Exceptions (%d) %r", len(failures), failures)
        mL = list(set(idList) - set(sL))
        logger.info("Missing (%d) %r", len(mL), mL)
        #
        return chainD, rawD

    def __imgtRemarkParser(self, pdbId, ifh):
        """IMGT REMARK 410 Parser

        Args:
            pdbId (str): input PDB ID [description]
            ifh (obj): input file handle

        Returns:
            dict: content dictionary of parsed details
        """
        sD = {
            "IMGT protein name": {"section": "proteins"},
            "ligand(s)": {"section": "ligands"},
            "Chain ID  ": {"section": "chains"},
        }
        pD = {
            "Chain ID  ": {"ky": "chain_data", "action": "appendAll"},
            #
            "ligand(s)": {"ky": "ligands", "action": "appendLine"},
            "IMGT protein name": {"ky": "proteinName", "action": "appendLine"},
            "IMGT receptor type": {"ky": "receptorType", "action": "appendLine"},
            "IMGT receptor description": {"ky": "receptorDescription", "action": "appendLine"},
            "Species": {"ky": "species", "action": "appendLine"},
            "Chain ID": {"ky": "chain_ids", "action": "appendLine"},
            #
        }
        cD = {}
        oD = {}
        curSection = None
        action = None
        curKy = None
        curSection = None
        curChain = None
        for ul in ifh.readlines():
            line = ul.decode("utf-8")
            if not line.startswith("REMARK 410 "):
                continue
            #
            curLine = line[11:-1]

            for section, sectionD in sD.items():
                if curLine.startswith(section):
                    logger.debug("%r Detected section %r", pdbId, section)
                    curSection = sectionD["section"]
                    first = True
                    break
            #
            for label, labelD in pD.items():
                if curLine.startswith(label):
                    logger.debug("%r detected label %r", pdbId, label)
                    curKy = labelD["ky"]
                    action = labelD["action"]
                    first = True
                    break
            #
            if action == "appendLine":
                if first:
                    first = False
                    logger.debug("Skipped %r", curKy)
                    continue
                logger.debug(">> SECTION %r KEY %r Adding %r", curSection, curKy, curLine.strip())
                oD.setdefault(curSection, {}).setdefault(curKy, []).append(curLine.strip())
            elif action == "appendAll":
                if first:
                    tL = [t for t in curLine.split(" ") if t]
                    curChain = tL[2]
                    first = False
                    logger.debug("%r current chain key %r", pdbId, curChain)
                    continue
                oD.setdefault(curSection, {}).setdefault(curChain, []).append(curLine)
        # --  rD raw extracted REM 410 content
        #  Post-process the domain annotations and alignments

        for chId, cL in oD["chains"].items() if "chains" in oD else {}:
            logger.debug("%r chainId %r (%d)", pdbId, chId, len(cL))
            tD = {}
            tD["description"] = self.__getField(cL, label="IMGT chain description  ")
            tD["domains"] = self.__splitDomains(pdbId, cL)
            #
            aD = self.__getAlignment(pdbId, cL)
            if aD and (len(aD["alignMapDL"]) == len(tD["domains"])):
                aL = aD["alignMapDL"]
                for ii, dD in enumerate(tD["domains"].values()):
                    if aL and len(aL) > ii:
                        dD["alignment"] = aL[ii]
            # --
            #  Integrate raw "proteins"  content
            if "proteins" in oD:
                paD = self.__getProteinAnnotations(chId, oD["proteins"])
                logger.debug("paD %r", paD)
                tD.update(paD)
            # --
            cD[chId] = tD
        #
        return cD, oD
        #

    def __getProteinAnnotations(self, chainId, pLD):
        """
        Example:

         "proteins": {
            "proteinName": [
               "IgG4 Sigma1 Fc"
            ],
            "receptorType": [
               "IG"
            ],
            "receptorDescription": [
               "FUSION-[TNFRSF1B]2-FC-GAMMA-1"
            ],
            "species": [
               "Homo sapiens (human)"
            ],
            "chain_ids": [
               "5w5m_A,5w5m_B"
            ]
         },
        """
        retD = {}
        try:
            ind = -1
            if "chain_ids" in pLD:
                for ii, chS in enumerate(pLD["chain_ids"]):
                    if chainId in chS:
                        ind = ii
                        break
                if ind >= 0:
                    for ky in ["proteinName", "receptorType", "receptorDescription", "species"]:
                        if ky in pLD and len(pLD[ky]) > ind:
                            retD[ky] = pLD[ky][ind]
                else:
                    logger.info("missing chain %r in %r", chainId, pLD["chain_ids"])
            else:
                logger.info("missing chain details for %r in %r", chainId, pLD)
        except Exception as e:
            logger.exception("Failing for %r with %s", chainId, str(e))

        return retD

    def __getField(self, lineList, label):
        label = "IMGT chain description  "
        ret = None
        for line in lineList:
            if line.startswith(label):
                ret = line[len(label) :]
                break
        return ret

    def __splitDomains(self, pdbId, lineList):
        retD = {}
        startLabel1 = "-DOMAIN      IMGT domain description  "
        startLabel2 = "-LIKE-DOMAIN IMGT domain description  "
        startLabel3 = "-LIKE-DOMAIN IMGT domain description "
        #
        geneLabel1 = "-DOMAIN      IMGT gene and allele     "
        geneLabel2 = "-LIKE-DOMAIN IMGT gene and allele     "
        geneLabel3 = "-LIKE-DOMAIN IMGT gene and allele    "
        tD = {}
        domain = None
        numD = 0
        for line in lineList:
            if line[1:].startswith(startLabel1):
                numD += 1
                domain = line.split(" ")[0].strip()
                description = line[len(startLabel1) + 1 :].strip()
                continue
            if line[1:].startswith(startLabel2):
                numD += 1
                domain = line.split(" ")[0].strip()
                description = line[len(startLabel2) + 1 :].strip()
                continue
            if line[2:].startswith(startLabel3):
                numD += 1
                domain = line.split(" ")[0].strip()
                description = line[len(startLabel3) + 1 :].strip()
                continue
            if domain and line.startswith(domain):
                tD.setdefault((domain, description, numD), []).append(line)
        #
        qD = {}
        for (domain, description, numD), cL in tD.items():
            for line in cL:
                if line[1:].startswith(geneLabel1):
                    qD.setdefault((domain + "|" + description + "|" + str(numD)), []).append(line[len(geneLabel1) + 1 :])
                if line[1:].startswith(geneLabel2):
                    qD.setdefault((domain + "|" + description + "|" + str(numD)), []).append(line[len(geneLabel2) + 1 :])
                if line[2:].startswith(geneLabel3):
                    qD.setdefault((domain + "|" + description + "|" + str(numD)), []).append(line[len(geneLabel3) + 1 :])
        #
        # "Homo sapiens IGHG4*01 (96.4%), Homo sapiens IGHG4*03 (96.4%), Homo sapiens IGHG4*04 (96.4%)",
        for ky, cL in qD.items():
            logger.debug("cL %r", cL)
            tS = "".join(cL)
            tS = " ".join(tS.split())
            logger.debug("tS %r", tS)
            #  handle some missing commas in the raw data -
            tS = tS.replace(") ", "),")
            logger.debug("TAX> %r tS %r", pdbId, tS)
            gnSL = tS.split(",")
            gDL = []
            for gnS in gnSL:
                tL = gnS.strip().split()
                logger.debug("tL %r", tL)
                geneAllele = tL[-2]
                taxName = " ".join(tL[:-2])
                gDL.append({"taxName": taxName, "geneAllele": geneAllele})
            retD[ky] = {"geneAlleles": gDL}

        return retD

    def __getAlignment(self, pdbId, lineList):
        try:
            startPat = "Chain amino acid sequence"
            endPat1 = "-DOMAIN"
            endPat2 = "-LIKE-DOMAIN"
            aL = []
            keep = False
            for line in lineList:
                if line.startswith(startPat):
                    keep = True
                    continue
                if line[1:].startswith(endPat1):
                    break
                if line[1:].startswith(endPat2):
                    break
                if line[2:].startswith(endPat2):
                    break
                if keep:
                    aL.append(line)
            #
            sLine = "".join(aL[1::2])
            mLine = "".join(aL[0::2])
            # Lots of cases with (UNK) sequences where REM 410 format is corrupt --
            # if "(UNK)" in sLine:
            #    logger.error("%r unknown or modified residue in one-letter-code sequence %r", pdbId, sLine[:30] + "...")
            #    return {}
            ok, indD = self.__findMatchingGroups(mLine, startGroup="[", endGroup="]")
            if not ok:
                logger.error("%r determining alignment boundaries fails", pdbId)
                return {}
            pdbRangeL = []
            for iBeg, iEnd in indD.items():
                if iEnd - iBeg <= 3:
                    continue
                pdbRangeL.append({"begEntitySeqId": iBeg + 1, "endEntitySeqId": iEnd + 1})
            ok, indD = self.__findMatchingGroups(mLine, startGroup="(", endGroup=")")
            if not ok:
                logger.error("%r determining alignment boundaries fails", pdbId)
                return {}
            imgtRangeL = []
            try:
                for k, v in indD.items():
                    tS = mLine[k + 1 : v]
                    tL = tS.split("-")
                    iBeg = tL[0]
                    iEnd = tL[1]
                    imgtRangeL.append({"begIMGTSeqId": iBeg, "endIMGTSeqId": iEnd})
            except Exception as e:
                logger.error("%r parsing boundaries fails with %r for %r", pdbId, str(e), mLine[:30] + "...")
                return {}
            #
            alignMapDL = []
            for pdbD, imgtD in zip(pdbRangeL, imgtRangeL):
                dD = pdbD
                dD.update(imgtD)
                alignMapDL.append(dD)
            return {"mapping": mLine, "pdbSeq": sLine, "alignMapDL": alignMapDL}
        except Exception as e:
            logger.exception("Failing %r with %r with %s", pdbId, mLine, str(e))
        return {}

    def __findMatchingGroups(self, strIn, startGroup="[", endGroup="]"):
        retD = {}
        dStack = []
        ok = True
        try:
            for i, cS in enumerate(strIn):
                if cS == startGroup:
                    dStack.append(i)
                elif cS == endGroup:
                    if len(dStack) == 0:
                        logger.error("No matching closing group at position: %r", str(i))
                        ok = False
                    retD[dStack.pop()] = i

            if len(dStack) > 0:
                logger.error("No matching opening group at: %r", strIn(dStack.pop()))
                ok = False
        except Exception:
            pass
        return ok, retD
