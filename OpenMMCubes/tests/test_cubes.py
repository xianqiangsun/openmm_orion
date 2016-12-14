import unittest

from OpenMMCubes.cubes import OpenMMComplexSetup
from OpenMMCubes.ports import OpenMMSystemInput, OpenMMSystemOutput
from simtk import openmm, unit

from floe.test import CubeTestRunner

from openeye import oechem


class SleepyCubeTester(unittest.TestCase):
    """
    Example test of an individual cube
    """

    def setUp(self):
        self.cube = OpenMMComplexSetup("sleepy")
        self.cube.args.receptor = 'OpenMMCubes/tests/input/receptor.pdbfixer.pdb'
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def tearDown(self):
        self.runner.finalize()

    def test_success(self):
        # Read a molecule
        mol = oechem.OEMol()
        mol2_filename = 'OpenMMCubes/tests/input/ligand.tripos.mol2'
        ifs = oechem.oemolistream(mol2_filename)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % mol2_filename)
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        self.assertEqual(self.runner.outputs["success"].qsize(), 1)

        # Check the output can be reconstituted into an OpenMM system object
        system = self.runner.outputs["success"].get()
        assert (system.getNumParticles() > 0)

    def test_failure(self):
        pass

    def test_ports(self):
        system = openmm.System()

        output = OpenMMSystemOutput("")
        encoded = output.encode(system)

        intake = OpenMMSystemInput("")
        decoded = intake.decode(encoded)

        assert (decoded.getNumParticles() == 0)
