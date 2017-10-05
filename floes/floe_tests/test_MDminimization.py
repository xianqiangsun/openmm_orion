from unittest import TestCase
from OpenMMCubes import utils
from floes.openmm_MDminimize import job as floe
import pytest


class MinimizationTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    @pytest.mark.slow
    def test_minimization(self):

        complex_path = utils.get_data_filename('examples', 'data/pbace_lcat13a_complex.oeb.gz')

        run_args = [
            '--system', complex_path,
            '--ofs-data_out', 'success.oeb',
            '--fail-data_out', 'fail.oeb',
            '--steps', '30000',
        ]

        self.assertFalse(floe.run(args=run_args))

