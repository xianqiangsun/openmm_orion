import io, os, time, traceback, base64, shutil
from openeye import oechem
import numpy as np
from simtk import unit, openmm
from simtk.openmm import app
import netCDF4 as netcdf

from floe.api import OEMolComputeCube, parameter, MoleculeInputPort, BinaryMoleculeInputPort, BinaryOutputPort, OutputPort, ParallelOEMolComputeCube
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES

from LigPrepCubes.ports import CustomMoleculeInputPort, CustomMoleculeOutputPort
import YankCubes.utils as utils
from YankCubes.utils import get_data_filename

import yank
from yank.yamlbuild import *
import textwrap
import subprocess

hydration_yaml_template = """\
---
options:
  minimize: no
  timestep: %(timestep)f*femtoseconds
  nsteps_per_iteration: %(nsteps_per_iteration)d
  number_of_iterations: %(number_of_iterations)d
  temperature: %(temperature)f*kelvin
  pressure: %(pressure)f*atmosphere
  anisotropic_dispersion_correction: no
  verbose: yes

molecules:
  input_molecule:
    # Don't change input.mol2
    filepath: input.mol2
    antechamber:
      charge_method: null

solvents:
  tip3p:
    nonbonded_method: PME
    nonbonded_cutoff: 9*angstroms
    clearance: 8*angstroms
  gbsa:
    nonbonded_method: NoCutoff
    implicit_solvent: OBC2
  vacuum:
    nonbonded_method: NoCutoff

systems:
  hydration-tip3p:
    solute: input_molecule
    solvent1: tip3p
    solvent2: vacuum
    leap:
      parameters: [leaprc.gaff, leaprc.protein.ff14SB, leaprc.water.tip3p]
  hydration-gbsa:
    solute: input_molecule
    solvent1: gbsa
    solvent2: vacuum
    leap:
      parameters: [leaprc.gaff, leaprc.protein.ff14SB, leaprc.water.tip3p]

protocols:
  protocol-tip3p:
    solvent1:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.00]
  protocol-gbsa:
    solvent1:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.00]
        lambda_sterics:        [1.00, 0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.00]
        lambda_sterics:        [1.00, 0.00]

experiments:
  system: hydration-%(solvent)s
  protocol: protocol-%(solvent)s
"""

class YankHydrationCube(ParallelOEMolComputeCube):
    title = "YankHydrationCube"
    description = """
    Compute the hydration free energy of a small molecule with YANK.

    This cube uses the YANK alchemical free energy code to compute the
    transfer free energy of one or more small molecules from gas phase
    to TIP3P solvent.

    See http://getyank.org for more information about YANK.
    """
    classification = ["Alchemical free energy calculations"]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 3600}, # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    # These can override YAML parameters
    nsteps_per_iteration = parameter.IntegerParameter('nsteps_per_iteration', default=500,
                                     help_text="Number of steps per iteration")

    timestep = parameter.DecimalParameter('timestep', default=2.0,
                                     help_text="Timestep (fs)")

    simulation_time = parameter.DecimalParameter('simulation_time', default=0.100,
                                     help_text="Simulation time (ns/replica)")

    temperature = parameter.DecimalParameter('temperature', default=300.0,
                                     help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter('pressure', default=1.0,
                                 help_text="Pressure (atm)")

    solvent = parameter.StringParameter('solvent', default='gbsa',
                                 help_text="Solvent choice: one of ['gbsa', 'tip3p']")

    def begin(self):
        # TODO: Is there another idiom to use to check valid input?
        if self.args.solvent not in ['gbsa', 'tip3p']:
            raise Exception("solvent must be one of ['gbsa', 'tip3p']")

        # Make substitutions to YAML here.
        # TODO: Can we override YAML parameters without having to do string substitutions?
        options = {
            'timestep' : self.args.timestep,
            'nsteps_per_iteration' : self.args.nsteps_per_iteration,
            'number_of_iterations' : int(np.ceil(self.args.simulation_time * unit.nanoseconds / (self.args.nsteps_per_iteration * self.args.timestep * unit.femtoseconds))),
            'temperature' : self.args.temperature,
            'pressure' : self.args.pressure,
            'solvent' : self.args.solvent,
        }

        self.yaml = hydration_yaml_template % options

        # Compute kT
        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA # Boltzmann constant
        self.kT = kB * (self.args.temperature * unit.kelvin)

    def process(self, mol, port):
        kT_in_kcal_per_mole = self.kT.value_in_unit(unit.kilocalories_per_mole)

        try:
            # Check that molecule is charged.
            is_charged = False
            for atom in mol.GetAtoms():
                if atom.GetPartialCharge() != 0.0:
                    is_charged = True
            if not is_charged:
                raise Exception('Molecule %s has no charges; input molecules must be charged.' % mol.GetTitle())

            # Write the specified molecule out to a mol2 file without changing its name.
            mol2_filename = 'input.mol2'
            ofs = oechem.oemolostream(mol2_filename)
            oechem.OEWriteMol2File(ofs, mol)

            # Undo oechem fuckery with naming mol2 substructures `<0>`
            from YankCubes.utils import unfuck_oechem_mol2_file
            unfuck_oechem_mol2_file(mol2_filename)

            # Run YANK on the specified molecule.
            from yank.yamlbuild import YamlBuilder
            yaml_builder = YamlBuilder(self.yaml)
            yaml_builder.build_experiments()

            # Analyze the hydration free energy.
            from yank.analyze import estimate_free_energies
            (Deltaf_ij_solvent, dDeltaf_ij_solvent) = estimate_free_energies(netcdf.Dataset('output/experiments/solvent1.nc', 'r'))
            (Deltaf_ij_vacuum,  dDeltaf_ij_vacuum)  = estimate_free_energies(netcdf.Dataset('output/experiments/solvent2.nc', 'r'))
            DeltaG_hydration = Deltaf_ij_vacuum[0,-1] - Deltaf_ij_solvent[0,-1]
            dDeltaG_hydration = np.sqrt(Deltaf_ij_vacuum[0,-1]**2 + Deltaf_ij_solvent[0,-1]**2)

            # Add result to original molecule
            oechem.OESetSDData(mol, 'DeltaG_yank_hydration', str(DeltaG_hydration * kT_in_kcal_per_mole))
            oechem.OESetSDData(mol, 'dDeltaG_yank_hydration', str(dDeltaG_hydration * kT_in_kcal_per_mole))

            # Emit molecule to success port.
            self.success.emit(mol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

        # clean up
        filenames_to_delete = ['input.mol2', 'output']
        for filename in filenames_to_delete:
            if os.path.exists(filename):
                try:
                    os.remove(filename)
                except:
                    shutil.rmtree(filename)
