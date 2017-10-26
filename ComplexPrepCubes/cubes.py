from ComplexPrepCubes import utils
from OpenMMCubes import utils as pack_utils
from floe.api import (OEMolComputeCube, ParallelOEMolComputeCube, parameter, MoleculeInputPort)
from openeye import oechem
import traceback
from simtk import unit
from simtk.openmm import app
from oeommtools import utils as oeommutils
from oeommtools.packmol import oesolvate
import parmed


class HydrationCube(ParallelOEMolComputeCube):
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

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10.0,
        help_text="Padding around protein for solvent box (angstroms)")

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=50.0,
        help_text="Salt concentration (millimolar)")

    ref_structure = parameter.BooleanParameter(
        'ref_structure',
        default=True,
        help_text="If Checked/True the molecule before solvation is attached to the solvated one")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, system, port):

        try:
            # Solvate the system
            sol_system = utils.hydrate(system, self.opt)
            sol_system.SetTitle(system.GetTitle())
            # Attached the original system to the solvated one
            if self.opt['ref_structure']:
                sol_system.SetData(oechem.OEGetTag("RefStructure"), system)
            self.success.emit(sol_system)
        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            system.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(system)

        return


class SolvationCube(ParallelOEMolComputeCube):
    title = "Solvation Cube Packmol"
    version = "0.0.0"
    classification = [["Preparation", "OEChem"]]
    tags = ['OEChem', 'PackMol']
    description = """
    This cube solvate the molecular system

    Input:
    -------
    oechem.OEMCMol - Streamed-in of the molecular system

    Output:
    -------
    oechem.OEMCMol - Emits the solvated system
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    density = parameter.DecimalParameter(
        'density',
        default=1.0,
        help_text="Solution density in g/ml")

    padding_distance = parameter.DecimalParameter(
        'padding_distance',
        default=10.0,
        help_text="The padding distance between the solute and the box edge in A")

    distance_between_atoms = parameter.DecimalParameter(
        'distance_between_atoms',
        default=2.0,
        help_text="The minimum distance between atoms in A")

    solvents = parameter.StringParameter(
        'solvents',
        required=True,
        default='[H]O[H]',
        help_text='Select solvents. The solvents are specified as comma separated smiles strings'
                  'e.g. [H]O[H], C(Cl)(Cl)Cl, CS(=O)C')

    molar_fractions = parameter.StringParameter(
        'molar_fractions',
        default='1.0',
        help_text="Molar fractions of each solvent components. The molar fractions are specified"
                  "as comma separated molar fractions strings e.g. 0.5,0.2,0.3")

    geometry = parameter.StringParameter(
        'geometry',
        default='box',
        choices=['box', 'sphere'],
        help_text="Geometry selection: box or sphere. Sphere cannot be used as periodic system "
                  "along with MD simulation")
    
    close_solvent = parameter.BooleanParameter(
        'close_solvent',
        default=False,
        help_text="If Checked/True solvent molecules will be placed very close to the solute")

    salt = parameter.StringParameter(
        'salt',
        default='[Na+], [Cl-]',
        help_text='Salt type. The salt is specified as list of smiles strings. '
                  'Each smiles string is the salt component dissociated in the '
                  'solution e.g. Na+, Cl-')

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=0.0,
        help_text="Salt concentration in millimolar")

    neutralize_solute = parameter.BooleanParameter(
        'neutralize_solute',
        default=True,
        help_text='Neutralize the solute by adding Na+ and Cl- counter-ions based on'
                  'the solute formal charge')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, solute, port):

        try:
            opt = dict(self.opt)
            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(solute) if dp.GetTag() in
                        ["solvents", "molar_fractions", "density"]}
            if new_args:
                for k in new_args:
                    if k == 'molar_fractions':
                        continue
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(solute.GetTitle(), new_args))
                opt.update(new_args)

            # Solvate the system
            sol_system = oesolvate(solute, **opt)
            self.log.info("Solvated System atom number {}".format(sol_system.NumAtoms()))
            sol_system.SetTitle(solute.GetTitle())
            self.success.emit(sol_system)
        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            solute.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(solute)

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

    remove_explicit_solvent = parameter.BooleanParameter(
        'remove_explicit_solvent',
        default=False,
        description='If True/Checked removes water and ion molecules from the system')

    system_port = MoleculeInputPort("system_port")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.wait_on('system_port')
        self.count = 0
        self.check_system = False

    def process(self, mol, port):
        try:
            if port == 'system_port':

                # Remove from solution water and ions
                if self.opt['remove_explicit_solvent']:
                    mol = oeommutils.strip_water_ions(mol)

                self.system = mol
                self.check_system = True
                return

            if self.check_system:
                num_conf = 0

                try:
                    lig_id = mol.GetData("IDTag")
                    name = 'p' + self.system.GetTitle() + '_' + lig_id
                except:
                    name = 'p' + self.system.GetTitle() + '_l' + mol.GetTitle()[0:12] + '_' + str(self.count)

                for conf in mol.GetConfs():
                    conf_mol = oechem.OEMol(conf)
                    complx = self.system.CreateCopy()
                    oechem.OEAddMols(complx, conf_mol)
                    
                    # Split the complex in components
                    protein, ligand, water, excipients = oeommutils.split(complx)

                    # If the protein does not contain any atom emit a failure
                    if not protein.NumAtoms():  # Error: protein molecule is empty
                        oechem.OEThrow.Fatal("The protein molecule does not contains atoms")

                    # If the ligand does not contain any atom emit a failure
                    if not ligand.NumAtoms():  # Error: ligand molecule is empty
                        oechem.OEThrow.Fatal("The Ligand molecule does not contains atoms")

                    # Check if the ligand is inside the binding site. Cutoff distance 3A
                    if not oeommutils.check_shell(ligand, protein, 3):
                        oechem.OEThrow.Fatal("The ligand is probably outside the protein binding site")

                    # Removing possible clashes between the ligand and water or excipients
                    if water.NumAtoms():
                        water_del = oeommutils.delete_shell(ligand, water, 1.5, in_out='in')

                    if excipients.NumAtoms():
                        excipient_del = oeommutils.delete_shell(ligand, excipients, 1.5, in_out='in')

                    # Reassemble the complex
                    new_complex = protein.CreateCopy()
                    oechem.OEAddMols(new_complex, ligand)
                    if excipients.NumAtoms():
                        oechem.OEAddMols(new_complex, excipient_del)
                    if water.NumAtoms():
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
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
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

    ligand_res_name = parameter.StringParameter(
        'ligand_res_name',
        required=True,
        default='LIG',
        help_text='Ligand residue name')

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
            protein, ligand, water, excipients = oeommutils.split(mol, ligand_res_name=self.opt['ligand_res_name'])

            self.log.info("\nComplex name: {}\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
                          "Water atom numbers = {}\nExcipients atom numbers = {}".format(mol.GetTitle(),
                                                                                         protein.NumAtoms(),
                                                                                         ligand.NumAtoms(),
                                                                                         water.NumAtoms(),
                                                                                         excipients.NumAtoms()))

            # Unique prefix name used to output parametrization files
            self.opt['prefix_name'] = mol.GetTitle()

            oe_mol_list = []
            par_mol_list = []

            # Apply FF to the Protein
            if protein.NumAtoms():
                oe_mol_list.append(protein)
                protein_structure = utils.applyffProtein(protein, self.opt)
                par_mol_list.append(protein_structure)

            # Apply FF to the ligand
            if ligand.NumAtoms():
                oe_mol_list.append(ligand)
                ligand_structure = utils.applyffLigand(ligand, self.opt)
                par_mol_list.append(ligand_structure)

            # Apply FF to water molecules
            if water.NumAtoms():
                oe_mol_list.append(water)
                water_structure = utils.applyffWater(water, self.opt)
                par_mol_list.append(water_structure)

            # Apply FF to the excipients
            if excipients.NumAtoms():
                excipient_structure = utils.applyffExcipients(excipients, self.opt)
                par_mol_list.append(excipient_structure)

                # The excipient order is set equal to the order in related
                # parmed structure to avoid possible atom index mismatching
                excipients = oeommutils.openmmTop_to_oemol(excipient_structure.topology,
                                                           excipient_structure.positions,
                                                           verbose=False)
                oechem.OEPerceiveBondOrders(excipients)
                oe_mol_list.append(excipients)

            # Build the overall Parmed structure
            complex_structure = parmed.Structure()

            for struc in par_mol_list:
                complex_structure = complex_structure + struc

            complx = oe_mol_list[0].CreateCopy()
            num_atom_system = complx.NumAtoms()

            for idx in range(1, len(oe_mol_list)):
                oechem.OEAddMols(complx, oe_mol_list[idx])
                num_atom_system += oe_mol_list[idx].NumAtoms()

            if not num_atom_system == complex_structure.topology.getNumAtoms():
                oechem.OEThrow.Fatal("Parmed and OE topologies mismatch atom number error")

            complx.SetTitle(mol.GetTitle())

            # Set Parmed structure box_vectors
            is_periodic = True
            try:
                vec_data = pack_utils.PackageOEMol.getData(complx, tag='box_vectors')
                vec = pack_utils.PackageOEMol.decodePyObj(vec_data)
                complex_structure.box_vectors = vec
            except:
                is_periodic = False
                self.log.warn("System has been parametrize without periodic box vectors for vacuum simulation")

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
                    # thisRes.SetName("LIG")
                    thisRes.SetHetAtom(True)
                oechem.OEAtomSetResidue(at, thisRes)

            if packed_complex.GetMaxAtomIdx() != complex_structure.topology.getNumAtoms():
                raise ValueError("OEMol complex and Parmed structure mismatch atom numbers")

            # Check if it is possible to create the OpenMM System
            if is_periodic:
                complex_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                               nonbondedCutoff=10.0 * unit.angstroms,
                                               constraints=app.HBonds,
                                               removeCMMotion=False)
            else:
                complex_structure.createSystem(nonbondedMethod=app.NoCutoff,
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
