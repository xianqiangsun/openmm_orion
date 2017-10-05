from unittest import TestCase
from OpenMMCubes import utils
from floes.openmm_MDprod import job as floe
import pytest


class ProdTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    @pytest.mark.slow
    def test_prod(self):

        complex_path = utils.get_data_filename('examples', 'data/pP38_lp38a_2x_complex.oeb.gz')
        run_args = [
            '--system', complex_path,
            '--picosec', '10',
            '--temperature', '300',
            '--pressure', '1',
            '--trajectory_interval', '10',
            '--reporter_interval', '100',
            '--suffix', 'prod',
            '--ofs-data_out', 'success.oeb',
            '--fail-data_out', 'fail.oeb',
        ]

        self.assertFalse(floe.run(args=run_args))
