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
# pylint: disable=line-too-long
import logging
import os.path
import time

from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.unichem import unichem_client

# from chembl_webresource_client.settings import Settings

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.seqalign.MiscUtils import MiscUtils


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
        try:
            if addTaxonomy:
                miscU = MiscUtils()
                outDirPath = os.path.join(cachePath, "uniprot_id_mapping_selected")
                taxMapFileName = "uniprot_taxonomy.pic"
                taxD = miscU.getUniprotXref(13, outDirPath, taxMapFileName, fmt="pickle", useCache=True)
            #
            fD = mU.doImport(chemblTargetRawPath, fmt="fasta", commentStyle="default")
            #
            for seqId, sD in fD.items():
                chemblId = seqId.strip().split(" ")[0].strip()
                unpId = seqId[seqId.find("[") + 1 : seqId.find("]")]
                seq = sD["sequence"]
                if addTaxonomy:
                    taxId = taxD[unpId] if unpId in taxD else "-1"
                    seqId = unpId + "|" + chemblId + "|" + taxId
                else:
                    seqId = unpId + "|" + chemblId
                uD.setdefault(unpId, []).append(chemblId)
                oD[seqId] = {"sequence": seq}
            #
            ok1 = mU.doExport(chemblTargetPath, oD, fmt="fasta")
            ok2 = mU.doExport(chemblTargetUniprotMapPath, uD, fmt="json")
            logger.info("ChEMBL mapping path %s (%d)", chemblTargetUniprotMapPath, len(uD))
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

    def getDrugData(self, moleculeChEMBLIdList):
        """Get drug data for the input ChEMBL molecule list.

        Args:
            moleculeChEMBLIdList (list): list of ChEMBL molecule identifiers

        Returns:
          (dict):  dictionary  {ChEMBId: {drug data}}

        """
        oD = {}
        chunkSize = 50
        try:
            for ii in range(0, len(moleculeChEMBLIdList), chunkSize):
                drug = new_client.drug  # pylint: disable=no-member
                drug.set_format("json")
                mDL = drug.filter(molecule_chembl_id__in=moleculeChEMBLIdList[ii : ii + chunkSize])
                if mDL:
                    logger.info("mDL (%d)", len(mDL))
                    for mD in mDL:
                        oD.setdefault(mD["molecule_chembl_id"], []).append(mD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD

    def getMoleculeData(self, moleculeChEMBLIdList):
        """Get molecule data for the input ChEMBL molecule list.

        Args:
            moleculeChEMBLIdList (list): list of ChEMBL molecule identifiers

        Returns:
          (dict):  dictionary  {ChEMBId: {molecule data}}

        """
        oD = {}
        chunkSize = 50
        try:
            for ii in range(0, len(moleculeChEMBLIdList), chunkSize):
                drug = new_client.molecule  # pylint: disable=no-member
                drug.set_format("json")
                mDL = drug.filter(molecule_chembl_id__in=moleculeChEMBLIdList[ii : ii + chunkSize])
                if mDL:
                    logger.info("mDL (%d)", len(mDL))
                    for mD in mDL:
                        oD.setdefault(mD["molecule_chembl_id"], []).append(mD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD

    def getUniChemData(self, inchiKeyList):
        """Get UniChem data for the input InChiKey list.

        Args:
            InChIList (list): list of InChI key molecule identifiers

        Returns:
          (dict):  dictionary  {ChEMBId: {molecule data}}

        """
        mapD = {
            1: {"name": "chembl", "baseUrl": "https://www.ebi.ac.uk/chembl/", "entryUrl": "https://www.ebi.ac.uk/chembldb/compound/inspect/"},
            3: {"name": "pdb", "baseUrl": "http://www.ebi.ac.uk/pdbe/", "entryUrl": "http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/"},
            2: {"name": "drugbank", "baseUrl": "http://drugbank.ca/", "entryUrl": "http://www.drugbank.ca/drugs/"},
            5: {"name": "pubchem_dotf", "baseUrl": "http://pubchem.ncbi.nlm.nih.gov/sources/sources.cgi", "entryUrl": "http://pubchem.ncbi.nlm.nih.gov/substance/"},
            4: {"name": "gtopdb", "baseUrl": "http://www.guidetopharmacology.org", "entryUrl": "http://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId="},
            11: {"name": "ibm", "baseUrl": "http://www-935.ibm.com/services/us/gbs/bao/siip/nih/", "entryUrl": "http://www-935.ibm.com/services/us/gbs/bao/siip/nih/?sid="},
            6: {"name": "kegg_ligand", "baseUrl": "http://www.genome.jp/kegg/ligand.html", "entryUrl": "http://www.genome.jp/dbget-bin/www_bget?"},
            9: {"name": "zinc", "baseUrl": "http://zinc15.docking.org", "entryUrl": "http://zinc15.docking.org/substances/"},
            8: {"name": "nih_ncc", "baseUrl": "http://nihsmr.evotec.com/evotec/", "entryUrl": ""},
            10: {"name": "emolecules", "baseUrl": "https://www.emolecules.com/", "entryUrl": "https://www.emolecules.com/cgi-bin/more?vid="},
            12: {"name": "atlas", "baseUrl": "http://www.ebi.ac.uk/gxa/home", "entryUrl": "http://www.ebi.ac.uk/gxa/query?conditionQuery="},
            7: {"name": "chebi", "baseUrl": "http://www.ebi.ac.uk/chebi/downloadsForward.do", "entryUrl": "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A"},
            14: {
                "name": "fdasrs",
                "baseUrl": "http://fdasis.nlm.nih.gov/srs/srs.jsp",
                "entryUrl": "http://fdasis.nlm.nih.gov/srs/ProxyServlet?mergeData=true&objectHandle=DBMaint&APPLICATION_NAME=fdasrs&actionHandle=default&nextPage=jsp/srs/ResultScreen.jsp&TXTSUPERLISTID=",
            },
            15: {"name": "surechembl", "baseUrl": "https://www.surechembl.org/search/", "entryUrl": "https://www.surechembl.org/chemical/"},
            21: {"name": "pubchem_tpharma", "baseUrl": "http://www.thomson-pharma.com/", "entryUrl": "http://pubchem.ncbi.nlm.nih.gov/substance/"},
            22: {"name": "pubchem", "baseUrl": "http://pubchem.ncbi.nlm.nih.gov", "entryUrl": "http://pubchem.ncbi.nlm.nih.gov/compound/"},
            27: {"name": "recon", "baseUrl": "https://vmh.uni.lu", "entryUrl": "https://vmh.uni.lu/"},
            28: {"name": "molport", "baseUrl": "https://www.molport.com/shop/index", "entryUrl": "https://www.molport.com/shop/molecule-link/"},
            31: {
                "name": "bindingdb",
                "baseUrl": "https://www.bindingdb.org/bind/index.jsp",
                "entryUrl": "http://www.bindingdb.org/bind/chemsearch/marvin/MolStructure.jsp?monomerid=",
            },
            41: {"name": "swisslipids", "baseUrl": "http://www.swisslipids.org/", "entryUrl": "http://www.swisslipids.org/"},
            29: {"name": "nikkaji", "baseUrl": "http://jglobal.jst.go.jp/en/", "entryUrl": "http://jglobal.jst.go.jp/en/redirect?Nikkaji_No="},
            32: {"name": "comptox", "baseUrl": "https://comptox.epa.gov/dashboard/", "entryUrl": "https://comptox.epa.gov/dashboard/"},
            33: {"name": "lipidmaps", "baseUrl": "http://www.lipidmaps.org", "entryUrl": "http://www.lipidmaps.org/data/LMSDRecord.php?LMID="},
            35: {"name": "carotenoiddb", "baseUrl": "http://carotenoiddb.jp/index.html", "entryUrl": "http://carotenoiddb.jp/Entries/"},
            36: {"name": "metabolights", "baseUrl": "http://www.ebi.ac.uk/metabolights/", "entryUrl": "http://www.ebi.ac.uk/metabolights/"},
            37: {"name": "brenda", "baseUrl": "https://www.brenda-enzymes.org/index.php", "entryUrl": "https://www.brenda-enzymes.org/ligand.php?brenda_ligand_id="},
            17: {"name": "pharmgkb", "baseUrl": "https://www.pharmgkb.org", "entryUrl": "https://www.pharmgkb.org/drug/"},
            18: {"name": "hmdb", "baseUrl": "http://www.hmdb.ca", "entryUrl": "http://www.hmdb.ca/metabolites/"},
            24: {
                "name": "nmrshiftdb2",
                "baseUrl": "http://nmrshiftdb.nmr.uni-koeln.de/portal/media-type/html/user/anon/page/default.psml/js_pane/P-Home",
                "entryUrl": "http://nmrshiftdb.org/molecule/",
            },
            25: {"name": "lincs", "baseUrl": "http://www.lincsproject.org/", "entryUrl": "http://identifiers.org/lincs.smallmolecule/"},
            39: {"name": "chemicalbook", "baseUrl": "https://www.chemicalbook.com", "entryUrl": "https://www.chemicalbook.com/ChemicalProductProperty_EN_"},
            20: {"name": "selleck", "baseUrl": "http://www.selleckchem.com", "entryUrl": "http://www.selleckchem.com/products/"},
            23: {"name": "mcule", "baseUrl": "https://mcule.com", "entryUrl": "https://mcule.com/"},
            26: {"name": "actor", "baseUrl": "https://actor.epa.gov", "entryUrl": "http://actor.epa.gov/actor/chemical.xhtml?casrn="},
            34: {"name": "drugcentral", "baseUrl": "http://drugcentral.org", "entryUrl": "http://drugcentral.org/drugcard/"},
            38: {"name": "rhea", "baseUrl": "http://www.rhea-db.org", "entryUrl": "http://www.rhea-db.org/searchresults?q=CHEBI:"},
        }
        oD = {}
        try:
            for ky in inchiKeyList:
                unc = unichem_client  # pylint: disable=no-member
                # unc.set_format("json")
                uDL = unc.get(ky)
                if uDL:
                    qD = {}
                    for uD in uDL:
                        if "src_id" in uD and int(uD["src_id"]) in mapD:
                            qD[mapD[int(uD["src_id"])]["name"]] = uD["src_compound_id"]
                    if qD:
                        oD[ky] = qD

        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD
