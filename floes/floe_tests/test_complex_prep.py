from unittest import TestCase
from OpenMMCubes import utils
from floes.openmm_complex_prep import job as floe


class ComplexPrepTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    def test_complex_prep(self):

        protein_path = utils.get_data_filename('examples', 'data/Bace_protein.pdb')
        ligand_path = utils.get_data_filename('examples', 'data/lig_CAT13a.oeb.gz')

        run_args = [
            '--protein', protein_path,
            '--Ligands-data_in', ligand_path,
            '--max_conformers', '800'
            '--ofs-data_out', 'success.oeb',
            '--fail-data_out', 'fail.oeb',
        ]

        self.assertFalse(floe.run(args=run_args))
