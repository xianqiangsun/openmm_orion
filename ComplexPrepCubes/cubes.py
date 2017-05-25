from ComplexPrepCubes import utils
from OpenMMCubes import utils as pack_utils
from floe.api import OEMolComputeCube, ParallelOEMolComputeCube, parameter, OEMolIStreamCube, MoleculeInputPort
from openeye import oechem
from LigPrepCubes import ff_utils


class Reader(OEMolIStreamCube):
    title = "Reader Cube"
    version = "0.0.0"
    classification = [["Reader Cube", "OEChem", "Reader Cube"]]
    tags = ['OEChem']
    description = """
        A Reader Cube 
        Input:
        -------
        oechem.OEMCMol - Streamed-in of the bio-molecular system.

        Output:
        -------
        oechem.OEMCMol - Emits the biomolecular system
        """
    protein_suffix = parameter.StringParameter(
        'protein_suffix',
        default='PRT',
        help_text='The protein suffix name to be used to identify the protein')

    def begin(self):
        self.opt = vars(self.args)

    def __iter__(self):
        ifs = oechem.oemolistream(self.opt['data_in'])
        for mol in ifs.GetOEMols():
            mol.SetTitle(self.opt['protein_suffix'])
            yield mol


class Splitter(OEMolComputeCube):
    title = "Splitter Cube"
    version = "0.0.0"
    classification = [["Protein Preparation", "OEChem", "Split Molecule"]]
    tags = ['OEChem', 'OEBio']
    description = """
        A starting  bio-molecular system is read in and split in three
        components: the protein, the ligand and the excipients by using
        OEBio functions

        Input:
        -------
        oechem.OEMCMol - Streamed-in of the bio-molecular system.

        Output:
        -------
        oechem.OEMCMol - Emits the molecules: protein, ligand
                         and excipients
        """

    def begin(self):
        self.opt = vars(self.args)

    def process(self, mol, port):

        # Split the biomolecular system
        protein, ligand, water, excipients = utils.split(mol)

        # self.log.info('Protein  number of atoms: {}'.format(protein.NumAtoms()))
        # self.log.info('Water number of atoms: {}'.format(water.NumAtoms()))
        # self.log.info('Excipients number of atoms: {}'.format(excipients.NumAtoms()))

        system = protein.CreateCopy()

        oechem.OEAddMols(system, water)
        oechem.OEAddMols(system, excipients)

        system.SetTitle(mol.GetTitle())

        #utils.order_check(system, 'splitter.log')

        # If the protein does not contain any atom emit a failure
        if not protein.NumAtoms():
            self.failure.emit(protein)
        else:
            self.success.emit(system)


class LigChargeCube(ParallelOEMolComputeCube):
    title = "Ligand Charge Cube"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem', 'Quacpac']
    description = """
           This cube charges the Ligand by using the AM1BCC method

           Input:
           -------
           oechem.OEMCMol - Streamed-in of the ligand molecules
           
           Output:
           -------
           oechem.OEMCMol - Emits the charge ligands
           """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 3600}, # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

    max_conformers = parameter.IntegerParameter(
        'max_conformers',
        default=800,
        help_text="Max number of ligand conformers")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        # self.count = 0

    def process(self, ligand, port):

        charged_ligand = None

        if not oechem.OEHasPartialCharges(ligand):
            # Ligand sanitation
            charged_ligand = ff_utils.sanitize(ligand)
            # Charge the ligand
            charged_ligand = ff_utils.assignELF10charges(charged_ligand,
                                                         max_confs=self.opt['max_conformers'],
                                                         strictStereo=True)

        # If the ligand has been charged then transfer the computed
        # charges to the starting ligand
        if charged_ligand:
            map_charges = {at.GetIdx():at.GetPartialCharge() for at in charged_ligand.GetAtoms()}
            for at in ligand.GetAtoms():
                at.SetPartialCharge(map_charges[at.GetIdx()])

        self.success.emit(ligand)


class SolvationCube(OEMolComputeCube):
    title = "Solvate Cube"
    version = "0.0.0"
    classification = [["Complex Preparation", "OEChem", "Complex preparation"]]
    tags = ['OEChem', 'OpenMM', 'PDBFixer']
    description = """
           This cube Solvate the passed molecule system

           Input:
           -------
           oechem.OEMCMol - Streamed-in of the molecule system

           Output:
           -------
           oechem.OEMCMol - Emits the solvated system
           """

    pH = parameter.DecimalParameter(
        'pH',
        default=7.4,
        help_text="Solvent pH used to select appropriate protein protonation state.")

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
        # Solvate the system
        sol_system = utils.solvate(system, self.opt)
        sol_system.SetTitle(system.GetTitle())
        #utils.order_check(sol_system, 'solvation.log')
        self.success.emit(sol_system)


class ComplexPrep(OEMolComputeCube):
    title = "Complex Cube Preparation"
    version = "0.0.0"
    classification = [["Complex Preparation", "OEChem", "Complex preparation"]]
    tags = ['OEChem']
    description = """
        The Complex made of the protein, ligand, internal water and excipients
        is assembled 
        
        Input:
        -------
        oechem.OEMCMol - Streamed-in of the system molecules made of
                         the protein, the ligand, the water and the 
                         excipients 
                         
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

        if port == 'system_port':
            self.system = mol
            self.check_system = True
            return

        if self.check_system:
            num_conf = 0
            name = 'p' + self.system.GetTitle() + '_l' + str(self.count)
            for conf in mol.GetConfs():
                conf_mol = oechem.OEMol(conf)
                complx = self.system.CreateCopy()
                complx = utils.delete_shell(conf_mol, complx, 1.5, in_out='in')
                oechem.OEAddMols(complx, conf_mol)
                name_c = name
                if mol.GetMaxConfIdx() > 1:
                    name_c = name + '_c' + str(num_conf)
                complx.SetData(oechem.OEGetTag('IDTag'), name_c)
                complx.SetTitle(name_c)
                num_conf += 1
                self.success.emit(complx)
            self.count += 1


class ForceFieldPrep(OEMolComputeCube):
    title = "Force Field Preparation"
    version = "0.0.0"
    classification = [["Force Field Preparation", "OEChem", "Force Field preparation"]]
    tags = ['OEChem', 'OEBio', 'OpenMM']
    description = """
        The complex is parametrized by using the selected force fields

        Input:
        -------
        oechem.OEMCMol - Streamed-in of the solvated complex
      
        Output:
        -------
        oechem.OEMCMol - Emits force field parametrized complex
        """

    protein_forcefield = parameter.DataSetInputParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        help_text='Forcefield parameters for protein')

    solvent_forcefield = parameter.DataSetInputParameter(
        'solvent_forcefield',
        default='tip3p.xml',
        help_text='Forcefield parameters for solvent')

    ligand_forcefield = parameter.StringParameter(
        'ligand_forcefield',
        required=True,
        default='GAFF2',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Forcefield to parameterize the ligand')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.ProcessProtein = True
        self.ProcessWater = True
        self.ProcessExcipients = True

    def process(self, mol, port):

        protein, ligand, water, excipients = utils.split(mol)

        if self.ProcessProtein:
            self.protein_structure = utils.applyffProtein(protein, self.opt)
            self.ProcessProtein = False

        if self.ProcessWater:
            self.water_structure = utils.applyffWater(water, self.opt)
            self.ProcessWater = False

        if self.ProcessExcipients:
            if excipients.NumAtoms() > 0:
                self.excipient_structure = utils.applyffExcipients(excipients, self.opt)
            self.ProcessExcipients = False

        self.ligand_structure = utils.applyffLigand(ligand, self.opt)

        if excipients.NumAtoms() > 0:
            complex_structure = self.protein_structure + self.ligand_structure + \
                                self.excipient_structure + self.water_structure

        else:
            complex_structure = self.protein_structure + self.ligand_structure + \
                                self.water_structure

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

        # Delete data from the complex
        #complx.DeleteData(oechem.OEGetTag('box_vectors'))

        # Attach the Parmed structure to the complex
        packed_complex = pack_utils.PackageOEMol.pack(complx, complex_structure)

        # Attach the reference positions to the complex
        ref_positions = complex_structure.positions
        packedpos = pack_utils.PackageOEMol.encodePyObj(ref_positions)
        packed_complex.SetData(oechem.OEGetTag('OEMDDataRefPositions'), packedpos)

        # ofs = oechem.oemolostream(packed_complex.GetData('IDTag')+'.oeb')
        # oechem.OEWriteMolecule(ofs, packed_complex)

        #utils.order_check(packed_complex, 'ff.log')

        self.success.emit(packed_complex)