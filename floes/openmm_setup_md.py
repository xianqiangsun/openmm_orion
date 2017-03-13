from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import ChargeMCMol, SMIRFFParameterization, GAFFParameterization, FREDDocking

job = WorkFloe("SetupMD")

job.description = """
**Set up an OpenMM complex for simulation and run 1000 steps of MD**

This floe will generate a fully solvated system with TIP3P, where the ligad is parameterized
with the SMIRFF forcefield parameters and the protein parameterized with amber99sbild.
It will then generate the OpenMM system to be run for 1000 steps of MD and
store the System and State for easy restarting of the simulation.
"""

job.classification = [
    ["OpenEye", "OpenMM"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="docked ligands")

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein')
complex_setup.promote_parameter('molecule_forcefield', promoted_name='molecule_forcefield')
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_forcefield')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_forcefield')

md_sim = OpenMMSimulation('md_sim')
md_sim.promote_parameter('complex_mol', promoted_name='complex_mol')
md_sim.set_parameters(complex_mol='complex.oeb.gz')

ofs = OEMolOStreamCube('ofs')
ofs.set_parameters(data_out="simulation.oeb.gz")

job.add_cubes(ifs, complex_setup, md_sim, ofs)
ifs.success.connect(complex_setup.intake)
complex_setup.success.connect(md_sim.intake)
md_sim.success.connect(ofs.intake)


if __name__ == "__main__":
    job.run()
