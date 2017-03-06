import os
from unittest import TestCase
from floes.openmm_complex_setup import job as complex_setup_floe


class OpenMMComplexSetupTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    def test_simple_run(self):
        protein_path = os.path.join(os.path.dirname(__file__), '..', '..', 'examples', 'data', 'T4-protein.pdb')
        ligand_path = os.path.join(os.path.dirname(__file__), '..', '..', 'examples', 'data', 'toluene.pdb')
        # run_args should be a list like sys.argv
        run_args = [
            '--protein', protein_path,
            '--ligand', ligand_path
        ]
        self.assertFalse(complex_setup_floe.run(args=run_args))
