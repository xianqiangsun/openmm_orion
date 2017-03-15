import io, os, traceback, string, random, parmed
import subprocess
import openmoltools
import tempfile
from openeye import oechem, oedocking, oeomega
from floe.api import (
    parameter, ParallelOEMolComputeCube, OEMolComputeCube, SinkCube, MoleculeInputPort,
    StringParameter, MoleculeOutputPort
)
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES

from LigPrepCubes.ports import (
    CustomMoleculeInputPort, CustomMoleculeOutputPort)
from OpenMMCubes.ports import ( ParmEdStructureInput, ParmEdStructureOutput,
    OpenMMSystemOutput, OpenMMSystemInput )

from smarty.forcefield import ForceField
from smarty.forcefield_utils import create_system_from_molecule
from OpenMMCubes.utils import download_dataset_to_file, get_data_filename

from LigPrepCubes import ff_utils
import OpenMMCubes.utils as utils


def _generateRandomID(size=5, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


class ChargeMCMol(OEMolComputeCube):
    title = "ChargeMCMol"
    classification = [["Ligand Preparation"]]
    tags = ['OpenMM', 'IDTagging', 'AtomTyping']
    description = """
    Uses openmoltools to perform the following:
    (1) 'normalize_molecule': checks aromaticity, add explicit hydrogens and renaming by IUPAC
    (2) Generate conformers with OMEGA
    (3) Assigns partial charges with oequacpac.OEAssignPartialCharges
    """

    def process(self, mol, port):
        try:
            if not mol.GetTitle():
                idtag = _generateRandomID()
                self.log.warn('Mol title not found, setting to {}'.format(idtag))
            else:
                # Store the IDTag from the SMILES file.
                idtag = mol.GetTitle()

            #Generate the charged molecule, keeping the first conf.
            charged_mol = ff_utils.assignCharges(mol, max_confs=800, strictStereo=True,
                                      normalize=True, keep_confs=1)
            # Store the IUPAC name from normalize_molecule
            iupac = [ charged_mol.GetTitle().strip() ]
            # Pack as list incase of commas in IUPUC

            # Keep it as an SD Tag
            oechem.OESetSDData(charged_mol, 'IUPAC', str(iupac))
            # Reset the charged mol title to the original IDTag
            charged_mol.SetTitle(idtag)
            charged_mol.SetData(oechem.OEGetTag('IDTag'), idtag)
            oechem.OESetSDData(charged_mol, 'IDTag', idtag)

            self.success.emit(charged_mol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

class LigandParameterization(OEMolComputeCube):
    title = "LigandParameterization"
    classification = [["Ligand Preparation"]]
    tags = ['OpenMM', 'Forcefields', 'GAFF', 'SMIRNOFF']
    description = """Parameterize the ligand with the chosen forcefield,
    Attach the ParmEd Structure to the OEMol.
    """

    molecule_forcefield = parameter.StringParameter(
        'molecule_forcefield',
        required=True,
        default='GAFF2',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Forcefield to parameterize the molecule')

    def process(self, mol, port):
        try:
            pmd = ff_utils.GenMolStructure(mol, self.args.molecule_forcefield)
            molecule_structure = pmd.parameterize()
            self.log.info(str(molecule_structure))

            oechem.OESetSDData(mol, 'NumAtoms', str(mol.NumAtoms()))
            oechem.OESetSDData(mol, 'Structure', str(molecule_structure))
            oechem.OESetSDData(mol, 'FF', str(self.args.molecule_forcefield) )
            mol.SetData(oechem.OEGetTag('IDTag'), mol.GetTitle())
            packedmol = utils.PackageOEMol.pack(mol, molecule_structure)

            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)


class FREDDocking(OEMolComputeCube):
    title = "FREDDocking"
    description = "Dock OE molecules"
    classification = [["Ligand Preparation"]]
    tags = ['Docking', 'FRED']
    #tags = [tag for lists in classification for tag in lists]

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    receptor = parameter.DataSetInputParameter(
        'receptor',
        required=True,
        help_text='Receptor OEB File')

    def begin(self):
        receptor = oechem.OEGraphMol()
        self.args.receptor = download_dataset_to_file(self.args.receptor)
        if not oedocking.OEReadReceptorFile(receptor, str(self.args.receptor)):
            raise Exception("Unable to read receptor from {0}".format(self.args.receptor))

        #Initialize Docking
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
        #return dockedMol

    def end(self):
        pass
