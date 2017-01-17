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
        pdbfilename = self.cube.args.protein
        #pdbfilename = 'OpenMMCubes/tests/input/toluene.pdb'
        ifs = oechem.oemolistream(pdbfilename)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % pdbfilename)
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

if __name__ == "__main__":
        unittest.main()
