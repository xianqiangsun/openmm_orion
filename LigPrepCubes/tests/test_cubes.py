from OpenMMCubes.ports import ( ParmEdStructureInput, ParmEdStructureOutput,
    OpenMMSystemOutput, OpenMMSystemInput )

import unittest, parmed, base64, pickle
from LigPrepCubes.ports import (CustomMoleculeInputPort, CustomMoleculeOutputPort)
from LigPrepCubes.cubes import ChargeMCMol, GAFFParameterization, SMIRFFParameterization, FREDDocking
import OpenMMCubes.utils as utils
from simtk import openmm, unit
from floe.test import CubeTestRunner
from openeye import oechem, oedocking

class ChargeMCMolTester(unittest.TestCase):
    """
    Test the OMEGA cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = ChargeMCMol('charge')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        ligand = utils.get_data_filename('examples','data/JF6_1.ism')

        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(ligand)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % ligand)
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check outmol has at least 1 conformer
        outmol = self.runner.outputs["success"].get()
        self.assertGreaterEqual(outmol.GetMaxConfIdx(), 1)

        # Loop through atoms and make sure partial charges were set
        for iat, oat in zip(mol.GetAtoms(), outmol.GetAtoms()):
            self.assertNotEqual(iat.GetPartialCharge(), oat.GetPartialCharge)

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()

class SMIRFFTester(unittest.TestCase):
    """
    Test the SMIRFF Parameterization cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = SMIRFFParameterization('smirff')
        self.cube.args.molecule_forcefield = utils.get_data_filename('smirff99Frosst','smirff99Frosst.ffxml')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        self.cube.args.ligand = utils.get_data_filename('examples','data/JF6_1-chgdmc.oeb.gz')
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(self.cube.args.ligand)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % self.cube.args.ligand)
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Get the output molecule
        outmol = self.runner.outputs["success"].get()

        # Check for the attached ParmEd Structure and IDTag
        self.assertTrue(outmol.HasData(oechem.OEGetTag('IDTag')))
        self.assertTrue(outmol.HasData(oechem.OEGetTag('Structure')))

        # Check that the Structure is parameterized
        encoded_structure = outmol.GetData(oechem.OEGetTag('Structure'))
        self.assertIsInstance(encoded_structure, str)
        decoded_structure = base64.b64decode(encoded_structure)
        struct_dict = pickle.loads(decoded_structure)
        struct = parmed.structure.Structure()
        struct.__setstate__(struct_dict)
        self.assertEqual(len(struct.atoms),mol.NumAtoms())
        struct_string = [x.strip() for x in str(struct).split(';')]
        self.assertIn('parametrized>', struct_string)

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()

class GAFFTester(unittest.TestCase):
    """
    Test the GAFF Parameterization cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = GAFFParameterization('gaff')
        self.cube.args.molecule_forcefield = 'GAFF2'
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        self.cube.args.ligand = utils.get_data_filename('examples','data/JF6_1-chgdmc.oeb.gz')
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(self.cube.args.ligand)
        if not oechem.OEReadMolecule(ifs, mol):
            raise Exception('Cannot read molecule from %s' % self.cube.args.ligand)
        ifs.close()

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)

        # Get the output molecule
        outmol = self.runner.outputs["success"].get()

        # Assert that one molecule was emitted on the success port
        #self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        #self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check for the ParmEd Structure
        encoded_structure = outmol.GetData(oechem.OEGetTag('Structure'))
        self.assertIsInstance(encoded_structure, str)
        decoded_structure = base64.b64decode(encoded_structure)
        struct_dict = pickle.loads(decoded_structure)
        struct = parmed.structure.Structure()
        struct.__setstate__(struct_dict)
        self.assertEqual(len(struct.atoms),mol.NumAtoms())
        print(struct)
        struct_string = [x.strip() for x in str(struct).split(';')]
        self.assertIn('parametrized>', struct_string)

        gafftmpfiles = ['ligand.frcmod', 'ligand.inpcrd', 'ligand.prmtop', 'ligand.gaff.mol2', 'tleap_commands', 'leap.log']
        utils.cleanup(gafftmpfiles)

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()

class FREDTester(unittest.TestCase):
    """
    Test the OMEGA cube
    Example inputs from `openmm_orion/examples/data`
    """
    def setUp(self):
        self.cube = FREDDocking('fred')
        self.cube.args.receptor = utils.get_data_filename('examples','data/epox_hydrolase_receptor.oeb.gz')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # Read a molecule
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(utils.get_data_filename('examples', 'data/JF6_1-smirff.oeb.gz'))
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
