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
from openmoltools.openeye import get_charges
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
            charged_mol = get_charges(mol, max_confs=800, strictStereo=True,
                                      normalize=True, keep_confs=-1)
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

class SMIRFFParameterization(OEMolComputeCube):
    title = "SMIRFFParameterization"
    classification = [["Ligand Preparation"]]
    tags = ['OpenMM', 'SMIRFF', 'Forcefields']
    description = """Parameterize the ligand with the SMIRNOFF parameters, Attach the ParmEd Structure to the OEMol.
    """

    molecule_forcefield = parameter.DataSetInputParameter(
        'molecule_forcefield',
        default='smirff99Frosst.ffxml',
        help_text='Forcefield FFXML file for molecule')

    def begin(self):
        try:
            with open(self.args.molecule_forcefield) as ffxml:
                self.mol_ff = ForceField(ffxml)
        except:
            raise RuntimeError('Error opening {}'.format(self.args.molecule_forcefield))

    def process(self, mol, port):
        try:
            mol_top, mol_sys, mol_pos = create_system_from_molecule(self.mol_ff, mol)
            molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys, xyz=mol_pos)
            molecule_structure.residues[0].name = "LIG"

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

class GAFFParameterization(OEMolComputeCube):
    title = "Attach GAFF parameters to OE molecules"
    description = """
    Parameterize a molecule with GAFF or GAFF2 parameters via AmberTools.
    Attach the resulting ParmEd Structure to the OEMol.
    """

    # TO DO: add ambermini to manifest/requirements/install as appropriate

    classification = [["Testing", "Ligand Preparation"]]
    tags = [tag for lists in classification for tag in lists]

    #molecule_forcefield = parameter.DataSetInputParameter(
    molecule_forcefield = parameter.StringParameter('molecule_forcefield',
        default='GAFF',
        help_text="GAFF forcefield to use: 'GAFF' or 'GAFF2'. Default: 'GAFF'")

    def begin(self):
        #Make sure here that selected forcefield is GAFF or GAFF2
        ff = self.args.molecule_forcefield
        if not ff in ['GAFF', 'GAFF2']:
            raise RuntimeError('Selected forcefield %s is not GAFF or GAFF2' % ff)

        #TO DO: Is there anything else I should do here to die early if this is going to fail?

        #Try to check if tleap is going to fail
        tleapin = """source leaprc.%s
quit
""" % ff.lower()
        file_handle = open('tleap_commands', 'w')
        file_handle.writelines(tleapin)
        file_handle.close()

        tmp = subprocess.getoutput('tleap -f tleap_commands')
        elements = tmp.split('\n')
        for elem in elements:
            if 'Could not open file' in elem:
                raise(RuntimeError('Error encountered trying to load %s in tleap.') % ff)
            if 'command not found' in elem:
                raise(RuntimeError('Error: requires tleap.'))


    def process(self, mol, port):
        ff = self.args.molecule_forcefield
        try:
            # Check that molecule is charged.
            is_charged = False
            for atom in mol.GetAtoms():
                if atom.GetPartialCharge() != 0.0:
                    is_charged = True
            if not is_charged:
                raise Exception('Molecule %s has no charges; input molecules must be charged.' % mol.GetTitle())


            # Determine formal charge (antechamber needs as argument)
            chg = 0
            for atom in mol.GetAtoms():
                chg+=atom.GetFormalCharge()

            # Write out mol to a mol2 file to process via AmberTools
            mol2file = tempfile.NamedTemporaryFile(suffix='.mol2')
            mol2filename = mol2file.name
            with oechem.oemolostream(mol2filename) as ofs:
                res = oechem.OEWriteConstMolecule(ofs, mol)
                if res != oechem.OEWriteMolReturnCode_Success:
                    raise RuntimeError("Error writing molecule %s to mol2." % mol.GetTitle())

            # Run antechamber to type and parmchk for frcmod
            # requires openmoltools 0.7.5 or later, which should be conda-installable via omnia
            gaff_mol2_filename, frcmod_filename = openmoltools.amber.run_antechamber( 'ligand', mol2filename, gaff_version = ff.lower(), net_charge = chg)

            # Run tleap using specified forcefield
            prmtop, inpcrd = openmoltools.amber.run_tleap('ligand', gaff_mol2_filename, frcmod_filename, leaprc = 'leaprc.%s' % ff.lower() )

            # Load via ParmEd
            molecule_structure = parmed.amber.AmberParm( prmtop, inpcrd )
            molecule_structure.residues[0].name = "LIG"

            # Pack parameters back into OEMol
            oechem.OESetSDData(mol, 'NumAtoms', str(mol.NumAtoms()))
            oechem.OESetSDData(mol, 'Structure', str(molecule_structure))
            oechem.OESetSDData(mol, 'FF', ff.upper() )
            mol.SetData(oechem.OEGetTag('IDTag'), mol.GetTitle())
            packedmol = utils.PackageOEMol.pack(mol, molecule_structure)
            self.success.emit(packedmol)

            # Emit
            self.success.emit(packedmol)

            # close file
            mol2file.close()

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

class OEBSinkCube(SinkCube):
    """
    A custom sink cube that writes molecules to a oeb.gz
    """
    #classification = [["Testing", "Output"]]
    title = "OEBSinkCube"
    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')

    directory = parameter.StringParameter('directory',
                                         default='output',
                                         description='Directory name')
    suffix = parameter.StringParameter('suffix',
                                        required=True,
                                        description='suffix to append')
    def begin(self):
        if not os.path.exists(self.args.directory):
            os.makedirs(self.args.directory)

    def write(self, mol, port):
        if 'idtag' in mol.GetData().keys():
            idtag = mol.GetData(oechem.OEGetTag('idtag'))
        outfname = '{}/{}-{}.oeb.gz'.format(self.args.directory,
                                           idtag, self.args.suffix)
        with oechem.oemolostream(outfname) as ofs:
            res = oechem.OEWriteConstMolecule(ofs, mol)
            if res != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Error writing {}.oeb.gz".format(outfname))
            else:
                self.log.info('Saving to {}'.format(outfname))
