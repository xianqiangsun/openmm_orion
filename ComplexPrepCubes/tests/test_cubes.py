import unittest, parmed, base64, pickle
import ComplexPrepCubes.utils as utils
import OpenMMCubes.utils
from simtk import openmm, unit
from simtk.openmm import app
from floe.test import CubeTestRunner
from openeye import oechem


class ConversionTester(unittest.TestCase):
    """
    Test the Complex assemblecube
    Example inputs from `openmm_orion/examples/data/T4-protein.pdb`
    """
    def setUp(self):
        pass

    def test_oemol_to_openmmTop(self):
        protein = OpenMMCubes.utils.get_data_filename('examples','data/T4-protein.pdb')
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(protein)
        oechem.OEReadMolecule(ifs, mol)

        top, omm_pos = utils.oemol_to_openmmTop(mol)

        # Assert 
        self.assertEqual(top.getNumAtoms(), mol.NumAtoms())

        for (op_at, oe_at) in zip(top.atoms(), mol.GetAtoms()):
            self.assertEqual(op_at.index, oe_at.GetIdx())


        oe_pos = [v for k,v in mol.GetCoords().items()]
        self.assertEqual(oe_pos, omm_pos.in_units_of(unit.angstrom)/unit.angstrom)



    def test_openmmTop_to_oemol(self):
        protein = OpenMMCubes.utils.get_data_filename('examples','data/T4-protein.pdb')
        
        pdb = app.PDBFile(protein)

        oe_mol = utils.openmmTop_to_oemol(pdb.topology, pdb.positions)

        # Assert 
        self.assertEqual(pdb.topology.getNumAtoms(), oe_mol.NumAtoms())

        for (op_at, oe_at) in zip(pdb.topology.atoms(), oe_mol.GetAtoms()):
            self.assertEqual(op_at.index, oe_at.GetIdx())

        import numpy as np
        oe_pos = [v for k,v in oe_mol.GetCoords().items()]
        np.testing.assert_almost_equal(pdb.getPositions(asNumpy=True).in_units_of(unit.angstrom)/unit.angstrom,np.array(oe_pos),decimal=2)
     
        
if __name__ == "__main__":
        unittest.main()
