from ComplexPrepCubes import utils
from OpenMMCubes import utils as pack_utils
from floe.api import (OEMolComputeCube, ParallelOEMolComputeCube, parameter,
                      OEMolIStreamCube, MoleculeInputPort)
from openeye import oechem
import traceback
from simtk import unit
from simtk.openmm import app


class ProteinReader(OEMolIStreamCube):
    title = "Protein Reader Cube"
    version = "0.0.0"
    classification = [["Reader Cube", "OEChem", "Reader Cube"]]
    tags = ['OEChem']
    description = """
        A Protein Reader Cube 
        Input:
        -------
        oechem.OEMCMol or - Streamed-in of the bio-molecular system
        The input file can be an .oeb, .oeb.gz, .pdb  or a .mol2 file
        
        Output:
        -------
        oechem.OEMCMol - Emits the bio-molecular system
        """
    protein_prefix = parameter.StringParameter(
        'protein_prefix',
        default='PRT',
        help_text='The protein suffix name used to identify the protein')

    def begin(self):
        self.opt = vars(self.args)

    def __iter__(self):
        for mol in super(ProteinReader, self).__iter__():
            mol.SetTitle(self.args.protein_prefix)
            yield mol


class Splitter(OEMolComputeCube):
    title = "Splitter Cube"
    version = "0.0.0"
    classification = [["Protein Preparation", "OEChem", "Split Molecule"]]
    tags = ['OEChem', 'OEBio']
    description = """
        A starting  bio-molecular system is read in and split in four
        components: the protein, the ligand, the water and the excipients 
        by using OEBio functions

        Input:
        -------
        oechem.OEMCMol - Streamed-in of the bio-molecular system.

        Output:
        -------
        oechem.OEMCMol - Emits the assembled system made by the
         protein,  water and excipients
        """

    def begin(self):
        self.opt = vars(self.args)

    def process(self, mol, port):

        try:
            # Split the bio molecular system
            protein, ligand, water, excipients = utils.split(mol)

            # self.log.info('Protein  number of atoms: {}'.format(protein.NumAtoms()))
            # self.log.info('Ligand  number of atoms: {}'.format(ligand.NumAtoms()))
            # self.log.info('Water number of atoms: {}'.format(water.NumAtoms()))
            # self.log.info('Excipients number of atoms: {}'.format(excipients.NumAtoms()))

            system = protein.CreateCopy()

            oechem.OEAddMols(system, water)
            oechem.OEAddMols(system, excipients)

            system.SetTitle(mol.GetTitle())

            # If the protein does not contain any atom emit a failure
            if not protein.NumAtoms():   # Error: protein molecule is empty
                oechem.OEThrow.Fatal("Splitting: Protein molecule after system splitting is empty")
            else:
                self.success.emit(system)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)

        return


class SolvationCube(OEMolComputeCube):
    title = "Solvation Cube"
    version = "0.0.0"
    classification = [["Complex Preparation", "OEChem", "Complex preparation"]]
    tags = ['OEChem', 'OpenMM', 'PDBFixer']
    description = """
           This cube solvate the molecular system

           Input:
           -------
           oechem.OEMCMol - Streamed-in of the molecular system

           Output:
           -------
           oechem.OEMCMol - Emits the solvated system
           """

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10,
        help_text="Padding around protein for solvent box (angstroms)")

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=100,
        help_text="Salt concentration (millimolar)")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, system, port):

        try:
            # Solvate the system
            sol_system = utils.solvate(system, self.opt)
            sol_system.SetTitle(system.GetTitle())
            self.success.emit(sol_system)
        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            system.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(system)

        return


class ComplexPrep(OEMolComputeCube):
    title = "Complex Preparation Cube"
    version = "0.0.0"
    classification = [["Complex Preparation", "OEChem", "Complex preparation"]]
    tags = ['OEChem']
    description = """
        This cube assembles the complex made of the solvated system and the 
        ligands. If a ligand presents multiple conformers, then each conformer 
        is bonded to the protein to form the solvated complex. For example if a 
        ligand has 3 conformers then 3 complexes are generated.
        
        Input:
        -------
        oechem.OEMCMol - Streamed-in of the solvated system and the ligands
                         
        Output:
        -------
        oechem.OEMCMol - Emits the complex molecule
        """

    system_port = MoleculeInputPort("system_port")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.wait_on('system_port')
        self.count=0
        self.check_system = False

    def process(self, mol, port):
        try:
            if port == 'system_port':
                self.system = mol
                self.check_system = True
                return

            if self.check_system:
                num_conf = 0
                name = 'p' + self.system.GetTitle() + '_l' + mol.GetTitle()[0:12] + '_' + str(self.count)
                for conf in mol.GetConfs():
                    conf_mol = oechem.OEMol(conf)
                    complx = self.system.CreateCopy()
                    oechem.OEAddMols(complx, conf_mol)

                    # Split the complex in components
                    protein, ligand, water, excipients = utils.split(complx)

                    # If the protein does not contain any atom emit a failure
                    if not protein.NumAtoms():  # Error: protein molecule is empty
                        oechem.OEThrow.Fatal("The protein molecule does not contains atoms")

                    # If the ligand does not contain any atom emit a failure
                    if not ligand.NumAtoms():  # Error: ligand molecule is empty
                        oechem.OEThrow.Fatal("The Ligand molecule does not contains atoms")

                    # If the water does not contain any atom emit a failure
                    if not water.NumAtoms():  # Error: water molecule is empty
                        oechem.OEThrow.Fatal("The water does not contains atoms. This could happen if not"
                                             "solvation process has occurred")

                    # Check if the ligand is inside the binding site. Cutoff distance 3A
                    if not utils.check_shell(ligand, protein, 3):
                        oechem.OEThrow.Fatal("The ligand is probably outside the protein binding site")

                    # Removing possible clashes between the ligand and water or excipients
                    water_del = utils.delete_shell(ligand, water, 1.5, in_out='in')
                    excipient_del = utils.delete_shell(ligand, excipients, 1.5, in_out='in')

                    # Reassemble the complex
                    new_complex = protein.CreateCopy()
                    oechem.OEAddMols(new_complex, ligand)
                    oechem.OEAddMols(new_complex, excipient_del)
                    oechem.OEAddMols(new_complex, water_del)

                    name_c = name
                    if mol.GetMaxConfIdx() > 1:
                        name_c = name + '_c' + str(num_conf)
                    new_complex.SetData(oechem.OEGetTag('IDTag'), name_c)
                    new_complex.SetTitle(name_c)
                    num_conf += 1
                    self.success.emit(new_complex)
                self.count += 1

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)

        return


class ForceFieldPrep(ParallelOEMolComputeCube):
    title = "Force Field Preparation Cube"
    version = "0.0.0"
    classification = [["Force Field Preparation", "OEChem", "Force Field preparation"]]
    tags = ['OEChem', 'OEBio', 'OpenMM']
    description = """
        Each complex is parametrized by using the selected force fields

        Input:
        -------
        oechem.OEMCMol - Streamed-in of complexes
      
        Output:
        -------
        oechem.OEMCMol - Emits force field parametrized complexes
        """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 3600}, # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

    protein_forcefield = parameter.DataSetInputParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        help_text='Force field parameters for protein')

    solvent_forcefield = parameter.DataSetInputParameter(
        'solvent_forcefield',
        default='tip3p.xml',
        help_text='Force field parameters for solvent')

    ligand_forcefield = parameter.StringParameter(
        'ligand_forcefield',
        required=True,
        default='GAFF2',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Force field to parametrize the ligand')

    other_forcefield = parameter.StringParameter(
        'other_forcefield',
        required=True,
        default='GAFF2',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Force field used to parametrize other molecules not recognized by the protein force field')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, mol, port):
        try:
            # Split the complex in components in order to apply the FF
            protein, ligand, water, excipients = utils.split(mol)

            # Unique prefix name used to output parametrization files
            self.opt['prefix_name'] = mol.GetTitle()

            # Apply FF to the Protein
            protein_structure = utils.applyffProtein(protein, self.opt)

            # Apply FF to water molecules
            water_structure = utils.applyffWater(water, self.opt)

            # Apply FF to the excipients
            if excipients.NumAtoms() > 0:
                excipient_structure = utils.applyffExcipients(excipients, self.opt)

                # The excipient order is set equal to the order in related
                # parmed structure to avoid possible atom index mismatching
                excipients = utils.openmmTop_to_oemol(excipient_structure.topology,
                                                      excipient_structure.positions)

            # Apply FF to the ligand
            ligand_structure = utils.applyffLigand(ligand, self.opt)

            # Build the Parmed structure
            if excipients.NumAtoms() > 0:
                complex_structure = protein_structure + ligand_structure + \
                                    excipient_structure + water_structure
            else:
                complex_structure = protein_structure + ligand_structure + water_structure

            num_atom_system = protein.NumAtoms() + ligand.NumAtoms() + excipients.NumAtoms() + water.NumAtoms()

            if not num_atom_system == complex_structure.topology.getNumAtoms():
                oechem.OEThrow.Fatal("Parmed and OE topologies mismatch atom number error")

            # Assemble a new OEMol complex in a specific order
            # to match the defined Parmed structure complex
            complx = protein.CreateCopy()
            oechem.OEAddMols(complx, ligand)
            oechem.OEAddMols(complx, excipients)
            oechem.OEAddMols(complx, water)

            complx.SetTitle(mol.GetTitle())

            # Set Parmed structure box_vectors
            vec_data = pack_utils.PackageOEMol.getData(complx, tag='box_vectors')
            vec = pack_utils.PackageOEMol.decodePyObj(vec_data)
            complex_structure.box_vectors = vec

            # Attach the Parmed structure to the complex
            packed_complex = pack_utils.PackageOEMol.pack(complx, complex_structure)

            # Attach the reference positions to the complex
            ref_positions = complex_structure.positions
            packedpos = pack_utils.PackageOEMol.encodePyObj(ref_positions)
            packed_complex.SetData(oechem.OEGetTag('OEMDDataRefPositions'), packedpos)

            # Set atom serial numbers, Ligand name and HETATM flag
            # oechem.OEPerceiveResidues(packed_complex, oechem.OEPreserveResInfo_SerialNumber)
            for at in packed_complex.GetAtoms():
                thisRes = oechem.OEAtomGetResidue(at)
                thisRes.SetSerialNumber(at.GetIdx())
                if thisRes.GetName() == 'UNL':
                    thisRes.SetName("LIG")
                    thisRes.SetHetAtom(True)
                oechem.OEAtomSetResidue(at, thisRes)

            if packed_complex.GetMaxAtomIdx() != complex_structure.topology.getNumAtoms():
                raise ValueError("OEMol complex and Parmed structure mismatch atom numbers")

            # Check if it is possible to create the OpenMM System
            system = complex_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                                    nonbondedCutoff=10.0 * unit.angstroms,
                                                    constraints=app.HBonds,
                                                    removeCMMotion=False)

            self.success.emit(packed_complex)
        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)

        return
