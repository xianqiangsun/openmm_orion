from unittest import TestCase
from OpenMMCubes import utils
from floes.openmm_MDnvt import job as floe
import pytest


class NvtTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    @pytest.mark.slow
    def test_nvt(self):

        complex_path = utils.get_data_filename('examples', 'data/pP38_lp38a_2x_complex.oeb.gz')
        run_args = [
            '--complex', complex_path,
            '--picosec', '10',
            '--temperature', '300',
            '--restraints', 'noh (ligand or protein)',
            '--restraintWt', '2.0',
            '--trajectory_interval', '10',
            '--reporter_interval', '100',
            '--suffix', 'nvt',
            '--ofs-data_out', 'success.oeb',
            '--fail-data_out', 'fail.oeb',
        ]

        self.assertFalse(floe.run(args=run_args))
