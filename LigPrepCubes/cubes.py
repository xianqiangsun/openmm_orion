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

def _generateRandomID(size=5, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

class SetIDTagfromTitle(OEMolComputeCube):
    title = "SetIDTagfromTitle"
    description = """
    Assigns IDtag from the OEMol's title via mol.GetTitle() or gets assigned a random one
    This cube is a place holder, it will eventually accomplish two things:
    (1) Assigns a unique title to used for naming output files
    (2) Split the different confs/poses from a single molecule while assigning each one a unique id tag.
    """
    classification = [["Testing", "Ligand Preparation"]]
    tags = [tag for lists in classification for tag in lists]

    def process(self, mol, port):
        #Check for OEMol title for ID labeling
        if not mol.GetTitle():
            idtag = _generateRandomID()
            oechem.OEThrow.Warning('No title found, setting to {}'.format(idtag))
            mol.SetTitle(idtag)
        else:
            idtag = mol.GetTitle()
        try:
            # Set AtomTypes
            oechem.OETriposAtomNames(mol)
            oechem.OETriposAtomTypeNames(mol)
            mol.AddData(oechem.OEGetTag('idtag'), mol.GetTitle())
            self.success.emit(mol)
        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)


class SMIRFFParameterization(OEMolComputeCube):
    title = "Attach FFXML to OE molecules"
    description = """
    Parameterize the ligand with the smirff99Frosst.ffxml parameters,
    which is parsed with smarty. Attach the System to the OEMol.
    """
    classification = [["Testing", "Ligand Preparation"]]
    tags = [tag for lists in classification for tag in lists]

    molecule_forcefield = parameter.DataSetInputParameter(
        'molecule_forcefield',
        default='smirff99Frosst.ffxml',
        help_text='Forcefield FFXML file for molecule')

    def begin(self):
        try:
            dfname = get_data_filename(self.args.molecule_forcefield)
            with open(dfname) as ffxml:
                self.mol_ff = ForceField(ffxml)
        except:
            raise RuntimeError('Error opening {}'.format(dfname))

    def process(self, mol, port):
        # Create a copy incase of error
        init_mol = oechem.OEMol(mol)
        try:
            mol_top, mol_sys, mol_pos = create_system_from_molecule(self.mol_ff, mol)
            molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys, xyz=mol_pos)
            molecule_structure.residues[0].name = "MOL"

            # Encode System/Structure, Attach to mol
            sys_out = OpenMMSystemOutput('sys_put')
            struct_out = ParmEdStructureOutput('struct_out')
            mol.SetData(oechem.OEGetTag('system'), sys_out.encode(mol_sys))
            mol.SetData(oechem.OEGetTag('structure'), struct_out.encode(molecule_structure))
            self.success.emit(mol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            init_mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(init_mol)

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
            # TO DO: Check that molecule HAS charges here (usually not having charges is a sign of a mistake)

            # Write out mol to a mol2 file to process via AmberTools
            mol2file = tempfile.NamedTemporaryFile(suffix='.mol2')
            mol2filename = mol2file.name
            with oechem.oemolostream(mol2filename) as ofs:
                res = oechem.OEWriteConstMolecule(ofs, mol)
                if res != oechem.OEWriteMolReturnCode_Success:
                    raise RuntimeError("Error writing molecule %s to mol2." % mol.GetTitle())

            # Run antechamber to type and parmchk for frcmod
            # requires openmoltools 0.7.5 or later, which should be conda-installable via omnia
            gaff_mol2_filename, frcmod_filename = openmoltools.amber.run_antechamber( 'ligand', mol2filename, gaff_version = ff.lower())

            # Run tleap using specified forcefield
            prmtop, inpcrd = openmoltools.amber.run_tleap('ligand', gaff_mol2_filename, frcmod_filename, leaprc = 'leaprc.%s' % ff.lower() )

            # Load via ParmEd
            molecule_structure = parmed.amber.AmberParm( prmtop, inpcrd )
            molecule_structure.residues[0].name = "LIG"

            # Pack parameters back into OEMol
            # Encode System/Structure, Attach to mol -- old/long way (waiting for #28)
            packedmol = oechem.OEMol(mol)
            sys_out = OpenMMSystemOutput('sys_put')
            struct_out = ParmEdStructureOutput('struct_out')
            packedmol.SetData(oechem.OEGetTag('structure'), struct_out.encode(molecule_structure))
            #Eventually this will look more like:
            #packedmol = utils.PackageOEMol.pack(mol, molecule_structure)

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
