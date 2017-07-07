import unittest
from LigPrepCubes.cubes import FREDDocking, LigChargeCube
import OpenMMCubes.utils as utils
from floe.test import CubeTestRunner
from openeye import oechem


class LigChargeTester(unittest.TestCase):
    """
    Test ELF10 charge cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = LigChargeCube('elf10charge')
        self.cube.args.max_conforms = 800
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # File name of a charged ligand
        lig_fname = utils.get_data_filename('examples', 'data/lig_CAT13a_chg.oeb.gz')

        # Read OEMol molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(lig_fname)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % lig_fname)
        ifs.close()

        mol_copy = mol.CreateCopy()
        # Set the partial charge to zero
        for at in mol_copy.GetAtoms():
            at.SetPartialCharge(0.0)

        # Process the molecules
        self.cube.process(mol_copy, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check outmol
        outmol = self.runner.outputs["success"].get()

        # Loop through atoms and make sure partial charges were set
        for iat, oat in zip(mol.GetAtoms(), outmol.GetAtoms()):
            self.assertNotEqual(iat.GetPartialCharge(), oat.GetPartialCharge)

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


class FREDTester(unittest.TestCase):
    """
    Test the FRED docking cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = FREDDocking('fred')
        self.cube.args.receptor = utils.get_data_filename('examples','data/T4-receptor.oeb.gz')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(utils.get_data_filename('examples', 'data/TOL-smnf.oeb.gz'))
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Get the output molecule, check that it has score.
        outmol = self.runner.outputs["success"].get()
        self.assertTrue(oechem.OEHasSDData(outmol, 'Chemgauss4'))

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()

if __name__ == "__main__":
        unittest.main()
