import os
from unittest import TestCase
import OpenMMCubes.utils as utils
from floes.smiles_complex_setup import job as complex_setup_floe


class OpenMMComplexSetupTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    def test_simple_run(self):
        protein_path =  utils.get_data_filename('examples', 'data/epox_hydrolase_apo-protein.pdb')
        receptor_path = utils.get_data_filename('examples', 'data/epox_hydrolase_receptor.oeb.gz')
        ligand_path = utils.get_data_filename('examples','data/JF6_1-smirff.oeb.gz')
        # run_args should be a list like sys.argv
        run_args = [
            '--protein', protein_path,
            '--ligand', ligand_path,
            '--receptor', receptor_path
        ]
        self.assertTrue(complex_setup_floe.run(args=run_args))
