import unittest
from ComplexPrepCubes.cubes import ComplexPrep, ForceFieldPrep
from openeye import oechem
from OpenMMCubes.cubes import utils as ommutils
from floe.test import CubeTestRunner


class ComplexPrepTester(unittest.TestCase):
    """
    Test the Complex Preparation  cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = ComplexPrep('ComplexPrep')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # File name
        fn_protein = ommutils.get_data_filename('examples', 'data/Bace_solvated.oeb.gz')
        fn_ligand = ommutils.get_data_filename('examples', 'data/lig_CAT13a_chg.oeb.gz')

        # Read Protein molecule
        protein = oechem.OEMol()

        with oechem.oemolistream(fn_protein) as ifs:
            oechem.OEReadMolecule(ifs, protein)

        # Read Ligand molecule
        ligand = oechem.OEMol()

        with oechem.oemolistream(fn_ligand) as ifs:
            oechem.OEReadMolecule(ifs, ligand)

        # Process the molecules
        self.cube.process(protein, self.cube.intake.name)
        # Why do I have to manually set these on?
        self.cube.check_system = True
        self.cube.system = protein
        self.cube.process(ligand, self.cube.system_port)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        complex = self.runner.outputs["success"].get()

        self.assertEquals(complex.GetMaxAtomIdx(), 52312)

    def tearDown(self):
        self.runner.finalize()

    def test_failure(self):
        pass


class ForceFieldPrepTester(unittest.TestCase):
    """
      Test the Complex Preparation  cube
      Example inputs from `openmm_orion/examples/data`
      """

    def setUp(self):
        self.cube = ForceFieldPrep('ForceFieldPrep')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_excipient_successGaff2(self):
        print('Testing cube:', self.cube.name)
        # File name
        fn_complex = ommutils.get_data_filename('examples',
                                                'data/pbace_lcat13a_solvated_complex.oeb.gz')

        # Read Protein molecule
        complex= oechem.OEMol()

        with oechem.oemolistream(fn_complex) as ifs:
            oechem.OEReadMolecule(ifs, complex)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'GAFF2'
        self.cube.args.other_forcefield = 'GAFF2'

        # Process the molecules
        self.cube.process(complex, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        complex = self.runner.outputs["success"].get()

    def test_excipient_successSmirnoff(self):
        print('Testing cube:', self.cube.name)
        # File name
        fn_complex = ommutils.get_data_filename('examples',
                                                'data/pbace_lcat13a_solvated_complex.oeb.gz')

        # Read Protein molecule
        complex = oechem.OEMol()

        with oechem.oemolistream(fn_complex) as ifs:
            oechem.OEReadMolecule(ifs, complex)

        # Selecting ligand and excipient parametrization
        self.cube.args.ligand_forcefield = 'SMIRNOFF'
        self.cube.args.other_forcefield = 'SMIRNOFF'

        # Process the molecules
        self.cube.process(complex, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    def test_protein_non_std_residue(self):
        print('Testing cube:', self.cube.name)
        # File name
        fn_complex = ommutils.get_data_filename('examples',
                                                'data/pCDK2_l1h1q_solvated_complex.oeb.gz')

        # Read Protein molecule
        complex = oechem.OEMol()

        with oechem.oemolistream(fn_complex) as ifs:
            oechem.OEReadMolecule(ifs, complex)

        # Process the molecules
        self.cube.process(complex, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

    def tearDown(self):
        self.runner.finalize()

    def test_failure(self):
        pass


if __name__ == "__main__":
        unittest.main()
