import traceback
from openeye import oechem, oedocking
import OpenMMCubes.utils as utils
from LigPrepCubes import ff_utils
from floe.api import OEMolComputeCube, ParallelOEMolComputeCube, parameter


class LigChargeCube(ParallelOEMolComputeCube):
    title = "Ligand Charge Cube"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem', 'Quacpac']
    description = """
           This cube charges the Ligand by using the ELF10 charge method

           Input:
           -------
           oechem.OEMCMol - Streamed-in of the ligand molecules

           Output:
           -------
           oechem.OEMCMol - Emits the charged ligands
           """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    max_conformers = parameter.IntegerParameter(
        'max_conformers',
        default=800,
        help_text="Max number of ligand conformers")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, ligand, port):

        try:
            charged_ligand = None

            # Ligand sanitation
            ligand = ff_utils.sanitizeOEMolecule(ligand)

            if not oechem.OEHasPartialCharges(ligand):
                # Charge the ligand
                charged_ligand = ff_utils.assignELF10charges(ligand,
                                                             self.opt['max_conformers'], strictStereo=True)

            # If the ligand has been charged then transfer the computed
            # charges to the starting ligand
            if charged_ligand:
                map_charges = {at.GetIdx(): at.GetPartialCharge() for at in charged_ligand.GetAtoms()}
                for at in ligand.GetAtoms():
                    at.SetPartialCharge(map_charges[at.GetIdx()])

            self.success.emit(ligand)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            ligand.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(ligand)

        return


class FREDDocking(OEMolComputeCube):
    title = "FRED Docking"
    version = "0.0.1"
    classification = [ ["Ligand Preparation", "OEDock", "FRED"],
    ["Ligand Preparation", "OEDock", "ChemGauss4"]]
    tags = ['OEDock', 'FRED']
    description = """
    Dock molecules using the FRED docking engine against a prepared receptor file.
    Return the top scoring pose.

    Input:
    -------
    receptor - Requires a prepared receptor (oeb.gz) file of the protein to dock molecules against.
    oechem.OEMCMol - Expects a charged multi-conformer molecule on input port.

    Output:
    -------
    oechem.OEMol - Emits the top scoring pose of the molecule with attachments:
        - SDData Tags: { ChemGauss4 : pose score }
    """

    receptor = parameter.DataSetInputParameter(
        'receptor',
        required=True,
        help_text='Receptor OEB File')

    def begin(self):
        receptor = oechem.OEGraphMol()
        self.args.receptor = utils.download_dataset_to_file(self.args.receptor)
        if not oedocking.OEReadReceptorFile(receptor, str(self.args.receptor)):
            raise Exception("Unable to read receptor from {0}".format(self.args.receptor))

        # Initialize Docking
        dock_method = oedocking.OEDockMethod_Hybrid
        if not oedocking.OEReceptorHasBoundLigand(receptor):
            oechem.OEThrow.Warning("No bound ligand, switching OEDockMethod to ChemGauss4.")
            dock_method = oedocking.OEDockMethod_Chemgauss4
        dock_resolution = oedocking.OESearchResolution_Default
        self.sdtag = oedocking.OEDockMethodGetName(dock_method)
        self.dock = oedocking.OEDock(dock_method, dock_resolution)
        if not self.dock.Initialize(receptor):
            raise Exception("Unable to initialize Docking with {0}".format(self.args.receptor))

    def clean(self, mol):
        mol.DeleteData('CLASH')
        mol.DeleteData('CLASHTYPE')
        mol.GetActive().DeleteData('CLASH')
        mol.GetActive().DeleteData('CLASHTYPE')

    def process(self, mcmol, port):
        try:
            dockedMol = oechem.OEMol()
            res = self.dock.DockMultiConformerMolecule(dockedMol, mcmol)
            if res == oedocking.OEDockingReturnCode_Success:
                oedocking.OESetSDScore(dockedMol, self.dock, self.sdtag)
                self.dock.AnnotatePose(dockedMol)
                score = self.dock.ScoreLigand(dockedMol)
                self.log.info("{} {} score = {:.4f}".format(self.sdtag, dockedMol.GetTitle(), score))
                oechem.OESetSDData(dockedMol, self.sdtag, "{}".format(score))
                self.clean(dockedMol)
                self.success.emit(dockedMol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mcmol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mcmol)

    def end(self):
        pass
