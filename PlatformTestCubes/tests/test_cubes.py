import unittest

from PlatformTestCubes.cubes import PlatformTestCube

from floe.test import CubeTestRunner

from openeye.oechem import OEMol, OESmilesToMol, OEExactGraphMatch


class SleepyCubeTester(unittest.TestCase):
    """
    Example test of an individual cube
    """

    def setUp(self):
        self.cube = PlatformTestCube("sleepy")
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def tearDown(self):
        self.runner.finalize()

    def test_success(self):
        mol = OEMol()
        OESmilesToMol(mol, str("c1(c(cccc1)OC(=O)C)C(=O)O"))
        self.cube.process(mol, self.cube.intake.name)
        self.assertEqual(self.runner.outputs["success"].qsize(), 1)
        self.assertEqual(self.runner.outputs["failure"].qsize(), 1)

        new_mol = self.runner.outputs["success"].get()
        self.assertTrue(OEExactGraphMatch(new_mol, mol))

    def test_failure(self):
        pass
