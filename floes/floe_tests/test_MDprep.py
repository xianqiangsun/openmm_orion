from unittest import TestCase
from OpenMMCubes import utils
from floes.openmm_MDprep import job as floe
import pytest


class MDPrepTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    @pytest.mark.slow
    def test_mdprep(self):

        protein_path = utils.get_data_filename('examples', 'data/Bace_protein.pdb')
        ligand_path = utils.get_data_filename('examples', 'data/lig_CAT13a.oeb.gz')

        run_args = [
            '--protein', protein_path,
            '--ligands', ligand_path,
            '--max_conformers', '800',
            '--protein_prefix', 'Bace',
            '--solvent_padding', '10',
            '--salt_conc', '100',
            '--protein_ff', 'amber99sbildn.xml',
            '--solvent_ff', 'tip3p.xml',
            '--ligand_ff', 'GAFF2',
            '--other_ff', 'GAFF2',

            '--min_steps', '30000',
            '--m_restraints', 'noh (ligand or protein)',
            '--m_restraintWt', '5.0',

            '--warm_psec', '10.0',
            '--w_restraints', 'noh (ligand or protein)',
            '--w_restraintWt', '2.0',
            '--w_trajectory_interval', '1000',
            '--w_reporter_interval', '10000',
            '--w_outfname', 'warmup',

            '--eq1_psec', '10.0',
            '--eq1_restraints', 'noh (ligand or protein)',
            '--eq1_restraintWt', '2.0',
            '--eq1_trajectory_interval', '1000',
            '--eq1_reporter_interval', '10000',
            '--eq1_outfname', 'equil1',

            '--eq2_psec', '10.0',
            '--eq2_restraints', 'noh (ligand or protein)',
            '--eq2_restraintWt', '0.5',
            '--eq2_trajectory_interval', '1000',
            '--eq2_reporter_interval', '10000',
            '--eq2_outfname', 'equil2',

            '--eq3_psec', '10.0',
            '--eq3_restraints', 'ca_protein or (noh ligand)',
            '--eq3_restraintWt', '0.1',
            '--eq3_trajectory_interval', '1000',
            '--eq3_reporter_interval', '10000',
            '--eq3_outfname', 'equil3',

            '--ofs-data_out', 'success.oeb',
            '--fail-data_out', 'fail.oeb',
        ]

        self.assertFalse(floe.run(args=run_args))
