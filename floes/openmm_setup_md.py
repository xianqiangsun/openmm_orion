from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation

job = WorkFloe("SetupMD")

job.description = """
**Set up OpenMM complex for simulation**

Check out the awesome stuff at the [OpenMM website](http://openmm.org)
"""

job.classification = [
    ["OpenEye", "OpenMM"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
# Promotes the parameter to something we can specify from the command line as "--ifs=..."
# Temporarily have ifs parameter overlap with complex_setup parameter
ifs.promote_parameter("data_in", promoted_name="ligand", description="docked ligands")

# the name of the object has to match the string: this is the name of myself
complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein')
complex_setup.promote_parameter('ligand', promoted_name='ligand')
complex_setup.promote_parameter('molecule_forcefield', promoted_name='molecule_forcefield')
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_forcefield')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_forcefield')

md_sim = OpenMMSimulation('md_sim')
md_sim.promote_parameter('state', promoted_name='state')
md_sim.promote_parameter('system', promoted_name='system')
md_sim.promote_parameter('complex_pdb', promoted_name='complex_pdb')

pdbin = OEMolIStreamCube("pdbin")
pdbin.set_parameters(data_in='combined_structure.pdb')

# this is hardwiring the filename to the molecules coming out of ofs
# note: requires decompression with lzma.decompress or gunzip system.xml.xz
system_save = FileOutputCube('system_save')
system_save.set_parameters(name="system.xml.xz")

state_save = FileOutputCube('state_save')
state_save.set_parameters(name="state.xml.xz")

ofs = FileOutputCube('ofs')
ofs.set_parameters(name='simulation.log')

job.add_cubes(ifs, complex_setup, system_save, pdbin, md_sim, state_save, ofs)
ifs.success.connect(complex_setup.intake)
complex_setup.success.connect(system_save.intake)


pdbin.success.connect(md_sim.intake)
md_sim.success.connect(ofs.intake)
md_sim.checkpoint.connect(state_save.intake)

if __name__ == "__main__":
    job.run()
