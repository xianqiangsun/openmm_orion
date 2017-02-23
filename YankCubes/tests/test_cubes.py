import unittest, os, parmed
from YankCubes.utils import get_data_filename
from YankCubes.cubes import YankHydrationCube
from simtk import openmm, unit
from floe.test import CubeTestRunner
from openeye import oechem


class YankHydrationCubeTester(unittest.TestCase):
    """
    Test the Yank hydration free energy calculation cube
    """
    def setUp(self):
        self.cube = YankHydrationCube("yank_hydration")
        self.cube.args.simulation_time = 0.005 # ns/replica
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        """
        Test that YankHydration can successfully process a single molecule.
        """
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(get_data_filename('freesolv_mini.oeb.gz'))

        # Test a single molecule
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        outmol = self.runner.outputs["success"].get()
        # Check that the number of atoms in input and output molecules match.
        self.assertEqual(outmol.NumAtoms(), mol.NumAtoms())
        # Check that a free energy of hydration has been attached
        self.assertTrue(oechem.OEHasSDData(outmol, 'DeltaG_yank_hydration'))
        self.assertTrue(oechem.OEHasSDData(outmol, 'dDeltaG_yank_hydration'))

    def test_failure(self):
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(get_data_filename('uranium-hexafluoride.sdf'))

        # Test a single molecule
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule')
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 0)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 1)

        outmol = self.runner.outputs["failure"].get()
        # Check that the number of atoms in input and output molecules match.
        self.assertEqual(outmol.NumAtoms(), mol.NumAtoms())
        # Check that an error message has been attached
        self.assertTrue(outmol.HasData(oechem.OEGetTag('error')))

    def test_ports(self):
        pass

    def tearDown(self):
        self.runner.finalize()


if __name__ == "__main__":
        unittest.main()
