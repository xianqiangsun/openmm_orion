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

hydration_yaml_default = """\
---
options:
  minimize: no
  verbose: no
  timestep: %(timestep)f*femtoseconds
  nsteps_per_iteration: %(nsteps_per_iteration)d
  number_of_iterations: %(number_of_iterations)d
  temperature: %(temperature)f*kelvin
  pressure: %(pressure)f*atmosphere

molecules:
  input_molecule:
    # Don't change input.mol2
    filepath: input.mol2
    # TODO: Can we autodetect whether molecule has charges or not?
    openeye:
      quacpac: am1-bcc
    antechamber:
      charge_method: null

solvents:
  pme:
    nonbonded_method: PME
    nonbonded_cutoff: 9*angstroms
    clearance: 16*angstroms
  GBSA:
    nonbonded_method: NoCutoff
    implicit_solvent: OBC2
  vacuum:
    nonbonded_method: NoCutoff

systems:
  hydration:
    solute: input_molecule
    solvent1: GBSA
    solvent2: vacuum
    leap:
      parameters: [leaprc.gaff, leaprc.protein.ff14SB, leaprc.water.tip3p]

protocols:
  hydration-protocol-explicit:
    solvent1:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.75, 0.50, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
        lambda_sterics:        [1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00]

  hydration-protocol-implicit:
    solvent1:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.00]
        lambda_sterics:        [1.00, 0.00]
    solvent2:
      alchemical_path:
        lambda_electrostatics: [1.00, 0.00]
        lambda_sterics:        [1.00, 0.00]

experiments:
  system: hydration
  protocol: hydration-protocol-implicit
"""

def run_cli(arguments):
    """Generic helper to run command line arguments"""
    # cli.main(argv=arguments.split())
    command = 'yank ' + arguments
    [stoutdata, sterrdata] = subprocess.Popen(command.split()).communicate()

    # TODO: Interpret suprocess data better
    if sterrdata:
        message = "An error return value (%s) was obtained:\n" % str(sterrdata)
        message += "\n"
        message += stoutdata
        message += "\n"
        raise Exception(message)

class YankHydrationCube(OEMolComputeCube):
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

    #Define Custom Ports to handle oeb.gz files
    intake = CustomMoleculeInputPort('intake')
    success = CustomMoleculeOutputPort('success')

    # These can override YAML parameters
    nsteps_per_iteration = parameter.IntegerParameter('nsteps_per_iteration', default=500,
                                     help_text="Number of steps per iteration")

    timestep = parameter.DecimalParameter('timestep', default=2.0,
                                     help_text="Timestep (fs)")

    simulation_time = parameter.DecimalParameter('simulation_time', default=0.005,
                                     help_text="Simulation time (ns/replica)")

    temperature = parameter.DecimalParameter('temperature', default=300.0,
                                     help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter('pressure', default=1.0,
                                 help_text="Pressure (atm)")

    # TODO: Check if this is the best way to present a large YAML file to be edited
    yaml_template = parameter.StringParameter('yaml',
                                        default=hydration_yaml_default,
                                        description='suffix to append')

    def begin(self):
        # Make substitutions to YAML here.
        # TODO: Can we override YAML parameters without having to do string substitutions?
        options = {
            'timestep' : self.args.timestep,
            'nsteps_per_iteration' : self.args.nsteps_per_iteration,
            'number_of_iterations' : int(np.ceil(self.args.simulation_time * unit.nanoseconds / (self.args.nsteps_per_iteration * self.args.timestep * unit.femtoseconds))),
            'temperature' : self.args.temperature,
            'pressure' : self.args.pressure,
        }
        self.yaml = self.args.yaml_template % options

        # Compute kT
        kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA # Boltzmann constant
        self.kT = kB * (self.args.temperature * unit.kelvin)

    def process(self, input_molecule, port):
        kT_in_kcal_per_mole = self.kT.value_in_unit(unit.kilocalories_per_mole)

        try:
            # Make a deep copy of the molecule to form the result molecule
            result_molecule = oechem.OEGraphMol(input_molecule)

            # Write the specified molecule out to a mol2 file
            # TODO: Can we read .oeb files directly into YANK?
            # TODO: Do we need to use a randomly-generated filename to avoid collisions?
            mol2_filename = 'input.mol2'
            ofs = oechem.oemolostream(mol2_filename)
            oechem.OEWriteMolecule(ofs, input_molecule)

            # Undo oechem fuckery with naming mol2 substructures `<0>`
            from YankCubes.utils import unfuck_oechem_mol2_file
            unfuck_oechem_mol2_file(mol2_filename)

            # Delete output directory if it already exists.
            if os.path.exists('output'):
                shutil.rmtree('output')

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
            result_molecule.SetData('DeltaG_hydration', DeltaG_hydration * kT_in_kcal_per_mole)
            result_molecule.SetData('dDeltaG_hydration', dDeltaG_hydration * kT_in_kcal_per_mole)
            self.success.emit(result_molecule)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            input_molecule.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(input_molecule)
