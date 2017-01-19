from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation

job = WorkFloe("SetupOpenMMComplex")

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


# this is hardwiring the filename to the molecules coming out of ofs
# note: requires decompression with lzma.decompress or gunzip system.xml.xz
system_save = FileOutputCube('system_save')
system_save.set_parameters(name="system.xml.xz")


job.add_cubes(ifs, complex_setup, system_save)
ifs.success.connect(complex_setup.intake)
complex_setup.success.connect(system_save.intake)

if __name__ == "__main__":
    job.run()
