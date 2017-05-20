from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMnvtCube

job = WorkFloe("NVT")

job.description = """
**NVT simulation of an OpenMM-ready solvated complex**

Ex. `data='examples/data'; python floes/openmm_MDnvt.py --complex $data/9PC1X-complex.oeb.gz --picosec 10`

Parameters:
-----------
complex (file): OEB file of the prepared protein:ligand complex

Optional:
--------
picosec (float): Number of picoseconds to warm up the complex.
temperature (decimal): target final temperature after warming.

Outputs:
--------
ofs: Outputs to a <idtag>-warmup.oeb.gz file
"""

job.classification = [['NVT']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")

nvt = OpenMMnvtCube('nvt')
nvt.promote_parameter('time', promoted_name='picosec')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(ifs, nvt, ofs, fail)
ifs.success.connect(nvt.intake)
nvt.success.connect(ofs.intake)
nvt.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
