import subprocess, tempfile, parmed
from openeye import oechem, oequacpac
import openmoltools
from openmoltools.openeye import *
from OpenMMCubes.utils import get_data_filename


def assignELF10charges(molecule, max_confs=800, strictStereo=True):
    """
     This function computes atomic partial charges for an OEMol by
     using the ELF10 method

    Parameters:
    -----------
    molecule : OEMol object
        The molecule that needs to be charged
    max_confs : integer
        The max number of conformers used to calculate the atomic partial charges
    strictStereo : bool
        a flag used to check if atoms need to have assigned stereo chemistry or not

    Return:
    -------
    mol_copy : OEMol
        a copy of the original molecule with assigned atomic partial charges
    """

    mol_copy = molecule.CreateCopy()

    # The passed molecule could have already conformers. If the conformer number
    # does not exceed the max_conf threshold then max_confs conformations will
    # be generated
    if not mol_copy.GetMaxConfIdx() > 200:
        # Generate up to max_confs conformers
        mol_copy = generate_conformers(mol_copy, max_confs=max_confs, strictStereo=strictStereo)

    # Assign MMFF Atom types
    if not oechem.OEMMFFAtomTypes(mol_copy):
        raise RuntimeError("MMFF atom type assignment returned errors")

    # ELF10 charges
    status = oequacpac.OEAssignCharges(mol_copy, oequacpac.OEAM1BCCELF10Charges())

    if not status:
        raise RuntimeError("OEAssignCharges returned error code %d" % status)

    return mol_copy


def sanitizeOEMolecule(molecule):
    """
    This function checks if the molecule has coordinates,
    explicit hydrogens and aromaticity. If the molecule
    does not have coordinates a fatal error is raised.
    If the molecule does not have hydrogens or aramatic
    flags are missing then a copy of the molecule is fixed

    Parameters:
    -----------
    molecule: OEMol
        The molecule to be checked

    Return:
    -------
    mol_copy: OEMol
        A copy of the checked molecule with fixed aromaticity
        and hydrogens
    """
    mol_copy = molecule.CreateCopy()

    # Check if the molecule has 3D coordinates
    if not oechem.OEGetDimensionFromCoords(mol_copy):
        oechem.OEThrow.Fatal("The molecule coordinates are set to zero")
    # Check if the molecule has hydrogens
    if not oechem.OEHasExplicitHydrogens(mol_copy):
        oechem.OEAddExplicitHydrogens(mol_copy)
    # Check if the molecule has assigned aromaticity
    if not mol_copy.HasPerceived(oechem.OEPerceived_Aromaticity):
        oechem.OEAssignAromaticFlags(mol_copy, oechem.OEAroModelOpenEye)

    # TEMPORARY PATCH FOR SMIRNOFF
    oechem.OETriposAtomNames(mol_copy)
    # # Check for any missing atom names, if found reassign all of them
    # if any([atom.GetName() == '' for atom in mol_copy.GetAtoms()]):
    #     oechem.OETriposAtomNames(mol_copy)

    return mol_copy


class ParamLigStructure(object):
    """
    Generates parametrized ParmEd structure of the molecule with a chosen force field

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The openeye molecule to be parameterized
    forcefield : str
        String specifying the forcefield parameters to be used
    prefix_name : str
        String specifying the output prefix filename

    Returns
    ---------
    packedmol : openeye.oechem.OEMol
        Openeye molecule with the ParmEd Structure attached.
    """

    def __init__(self, molecule, forcefield, prefix_name='ligand', delete_out_files=True):
        if not forcefield in ['SMIRNOFF', 'GAFF', 'GAFF2']:
            raise RuntimeError('Selected forcefield %s is not GAFF/GAFF2/SMIRNOFF' % forcefield)
        else:
            self.molecule = molecule
            self.forcefield = str(forcefield).strip()
            self.structure = None
            self.prefix_name = prefix_name
            self.delete_out_files = delete_out_files

    @staticmethod
    def checkTleap(self):
        # Try to check if tleap is going to fail
        with open('tleap_commands', 'w') as cmd:
            cmd.write( "source leaprc.%s; quit" % self.forcefield.lower() )
        tmp = subprocess.getoutput('tleap -f tleap_commands')
        elements = tmp.split('\n')
        for elem in elements:
            if 'Could not open file' in elem:
                raise RuntimeError('Error encountered trying to load %s in tleap.'% self.forcefield)
            if 'command not found' in elem:
                raise RuntimeError('Error: requires tleap.')
        return True

    def checkCharges(self, molecule):
        # Check that molecule is charged.
        is_charged = False
        for atom in molecule.GetAtoms():
            if atom.GetPartialCharge() != 0.0:
                is_charged = True
        if not is_charged:
            raise Exception('Molecule %s has no charges; input molecules must be charged.' % molecule.GetTitle())

    def getGaffStructure(self, molecule=None, forcefield=None):
        if not molecule:
            molecule = self.molecule
        if not forcefield:
            forcefield = self.forcefield

        #  Try to check if tleap is going to fail
        self.checkCharges(molecule)

        # Determine formal charge (antechamber needs as argument)
        chg = 0
        for atom in molecule.GetAtoms():
            chg+=atom.GetFormalCharge()

        # Write out mol to a mol2 file to process via AmberTools
        mol2file = tempfile.NamedTemporaryFile(suffix='.mol2')
        mol2filename = mol2file.name
        with oechem.oemolostream(mol2filename) as ofs:
            res = oechem.OEWriteConstMolecule(ofs, molecule)
            if res != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Error writing molecule %s to mol2." % molecule.GetTitle())

        # Run antechamber to type and parmchk for frcmod
        # requires openmoltools 0.7.5 or later, which should be conda-installable via omnia
        gaff_mol2_filename, frcmod_filename = openmoltools.amber.run_antechamber(self.prefix_name, mol2filename,
                                                                                 gaff_version=forcefield.lower(),
                                                                                 charge_method=None)

        # Run tleap using specified forcefield
        prmtop, inpcrd = openmoltools.amber.run_tleap(self.prefix_name, gaff_mol2_filename,
                                                      frcmod_filename,
                                                      leaprc='leaprc.%s' % forcefield.lower())

        # Load via ParmEd
        molecule_structure = parmed.amber.AmberParm(prmtop, inpcrd)

        if self.delete_out_files:
            os.remove(gaff_mol2_filename)
            os.remove(frcmod_filename)
            os.remove(prmtop)
            os.remove(inpcrd)

        return molecule_structure

    def parameterize(self):
        from openforcefield.utils.utils import generateSMIRNOFFStructure
        if self.forcefield == 'SMIRNOFF':
            structure = generateSMIRNOFFStructure(molecule)
        elif self.forcefield in ['GAFF', 'GAFF2']:
            structure = self.getGaffStructure()
        self.structure = structure
        return self.structure
