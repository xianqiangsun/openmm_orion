import unittest

from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from OpenMMCubes.ports import OpenMMSystemInput, OpenMMSystemOutput
from simtk import openmm, unit

from floe.test import CubeTestRunner

from openeye import oechem


class SetupCubeTester(unittest.TestCase):
    """
    Test the OpenMM complex setup cube
    """

    def setUp(self):
        self.cube = OpenMMComplexSetup("complex_setup")
        self.cube.args.protein = 'input/T4-protein.pdb'
        self.cube.args.ligand = 'input/toluene.pdb'
        self.cube.args.molecule_forcefield = 'input/forcefield/smirff99Frosst.ffxml'
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def tearDown(self):
        self.runner.finalize()

    def test_success(self):
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(self.cube.args.ligand)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % self.cube.args.ligand)
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



class SimulationCubeTester(unittest.TestCase):

    def setUp(self):
        self.cube = OpenMMSimulation("md_sim")
        self.runner = CubeTestRunner(self.cube)
        self.cube.args.complex_pdb = 'input/combined_structure.pdb'
        self.cube.args.system = 'input/system.xml.xz'
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def tearDown(self):
        self.runner.finalize()

    def test_success(self):
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(self.cube.args.complex_pdb)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % self.cube.args.complex_pdb)
        ifs.close()
        self.cube.process(mol, self.cube.intake.name)
        state = self.runner.outputs["checkpoint"].get()

        self.assertIs(type(state), openmm.State)

    def test_failure(self):
        pass

if __name__ == "__main__":
        unittest.main()
