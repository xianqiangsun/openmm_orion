import unittest
from ComplexPrepCubes.cubes import  SolvationCube, ComplexPrep, ForceFieldPrep
from simtk import unit
from simtk.openmm import app
from openeye import oechem
from oeommtools import utils as oeommtools
from ComplexPrepCubes import utils
from OpenMMCubes.cubes import utils as ommutils
from floe.test import CubeTestRunner
import numpy as np


class SolvationCubeTester(unittest.TestCase):
    """
    Test the Solvation cube
    Example inputs from `openmm_orion/examples/data`
    """

    def setUp(self):
        self.cube = SolvationCube('Solvation')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def test_success(self):
        print('Testing cube:', self.cube.name)
        # File name
        fname = ommutils.get_data_filename('examples', 'data/Bace_protein.pdb')

        # Read OEMol molecule
        mol = oechem.OEMol()

        with oechem.oemolistream(fname) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        # Process the molecules
        self.cube.process(mol, self.cube.intake.name)
        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        outmol = self.runner.outputs["success"].get()

        prot, lig, wat, other = utils.split(outmol)

        npa = prot.GetMaxAtomIdx()
        nla = lig.GetMaxAtomIdx()
        noa = other.GetMaxAtomIdx()
        nwa = wat.GetMaxAtomIdx()

        self.assertEquals(npa, 6044)
        self.assertEquals(nla, 0)
        # Ions added to excipients
        self.assertEquals(noa, 80)
        # Water molecules added
        self.assertEquals(nwa, 46152)

        # Check Box vectors
        box_vectors = outmol.GetData('box_vectors')
        box_vectors = ommutils.PackageOEMol.decodePyObj(box_vectors)
        box_vectors = box_vectors.in_units_of(unit.nanometers)

        self.assertAlmostEqual(box_vectors[0][0]/unit.nanometers, 8.23, delta=0.01)
        self.assertEqual(box_vectors[0][1]/unit.nanometers, 0.0)
        self.assertEqual(box_vectors[0][2] / unit.nanometers, 0.0)

        self.assertAlmostEqual(box_vectors[1][1]/unit.nanometers, 8.23, delta=0.01)
        self.assertEqual(box_vectors[1][0]/unit.nanometers, 0.0)
        self.assertEqual(box_vectors[1][2] / unit.nanometers, 0.0)

        self.assertAlmostEqual(box_vectors[2][2]/unit.nanometers, 8.23, delta=0.01)
        self.assertEqual(box_vectors[2][0]/unit.nanometers, 0.0)
        self.assertEqual(box_vectors[2][1] / unit.nanometers, 0.0)

    def tearDown(self):
        self.runner.finalize()

    def test_failure(self):
        pass


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

        #complex = self.runner.outputs["success"].get()

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


class ConversionTester(unittest.TestCase):
    """
    Test conversion functions between OE and OpenMM topologies
    Example inputs from `openmm_orion/examples/data/T4-protein.pdb`
    """
    def setUp(self):
        pass

    def test_oemol_to_openmmTop(self):
        protein = ommutils.get_data_filename('examples', 'data/T4-protein.pdb')
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(protein)
        oechem.OEReadMolecule(ifs, mol)

        top, omm_pos = oeommtools.oemol_to_openmmTop(mol)

        # Assert Atom numbers
        self.assertEqual(top.getNumAtoms(), mol.NumAtoms())

        for (op_at, oe_at) in zip(top.atoms(), mol.GetAtoms()):
            # Assert atom indexes
            self.assertEqual(op_at.index, oe_at.GetIdx())

        oe_pos = [v for k, v in mol.GetCoords().items()]
        # Assert atom positions
        self.assertEqual(oe_pos, omm_pos.in_units_of(unit.angstrom)/unit.angstrom)

        # Assert bond order
        dic_bond_openmm = {}
        for bond in top.bonds():
            # OpenMM atoms
            at0_idx = bond[0].index
            at1_idx = bond[1].index
            if at0_idx < at1_idx:
                dic_bond_openmm[(at0_idx, at1_idx)] = bond.order
            else:
                dic_bond_openmm[(at1_idx, at0_idx)] = bond.order

        dic_bond_oe = {}
        for bond in mol.GetBonds():
            # OE atoms
            at0_idx = bond.GetBgnIdx()
            at1_idx = bond.GetEndIdx()
            if at0_idx < at1_idx: 
                dic_bond_oe[(at0_idx, at1_idx)] = bond.GetOrder()
            else:
                dic_bond_oe[(at1_idx, at0_idx)] = bond.GetOrder()

        self.assertEqual(dic_bond_openmm, dic_bond_oe)

    def test_openmmTop_to_oemol(self):
        protein = ommutils.get_data_filename('examples', 'data/T4-protein.pdb')
        
        pdb = app.PDBFile(protein)

        oe_mol = oeommtools.openmmTop_to_oemol(pdb.topology, pdb.positions, verbose=False)

        # Assert 
        self.assertEqual(pdb.topology.getNumAtoms(), oe_mol.NumAtoms())

        for (op_at, oe_at) in zip(pdb.topology.atoms(), oe_mol.GetAtoms()):
            self.assertEqual(op_at.index, oe_at.GetIdx())

        oe_pos = [v for k, v in oe_mol.GetCoords().items()]
        np.testing.assert_almost_equal(pdb.getPositions(asNumpy=True).in_units_of(unit.angstrom)/unit.angstrom,
                                       np.array(oe_pos), decimal=2)
     
        
if __name__ == "__main__":
        unittest.main()
