from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, FileOutputCube
from OpenMMCubes.cubes import OpenMMComplexSetup

job = WorkFloe("SetupOpenMMComplex")

job.description = """
**Set up OpenMM complex for simulation**

Check out the awesome stuff at the [OpenMM website](http://openmm.org)
"""

job.classification = [
    ["OpenMM", "Protein-Ligand Complex Setup"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ofs = FileOutputCube("ofs")

# Promotes the parameter to something we can specify from the command line as "--ifs=..."
ifs.promote_parameter("data_in", promoted_name="ifs", description="posed ligand")

# this is hardwiring the filename to the molecules coming out of ofs
# note: requires decompression with lzma.decompress or gunzip system.xml.xz
ofs.set_parameters(name="system.xml.xz")

# the name of the object has to match the string: this is the name of myself
complex_setup = OpenMMComplexSetup("complex_setup")

job.add_cubes(ifs, ofs, complex_setup)
ifs.success.connect(complex_setup.intake)
complex_setup.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
