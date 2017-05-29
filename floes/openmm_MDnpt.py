from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMnptCube

job = WorkFloe("NPT")

job.description = """
**NPT simulation of an OpenMM-ready solvated complex**

Ex. `data='examples/data'; python floes/openmm_MDnpt.py --complex $data/9PC1X-complex.oeb.gz --picosec 10`

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

job.classification = [['NPT']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")

npt = OpenMMnptCube('npt')
npt.promote_parameter('time', promoted_name='picosec', default=10.0, description='Length of MD run in picoseconds')
npt.promote_parameter('restraints', promoted_name='restraints', default="ca_protein and (noh ligand)", description='Select mask to apply restarints')
npt.promote_parameter('restraintWt', promoted_name='restraintWt', default=2.0, description='Restraint weight')
npt.promote_parameter('trajectory_interval', promoted_name='trajectory_interval', default=10, description='Trajectory saving interval')
npt.promote_parameter('reporter_interval', promoted_name='reporter_interval', default=100, description='Reporter saving interval')
npt.promote_parameter('outfname', promoted_name='outfname', default='npt', description='Equilibration suffix name')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(ifs, npt, ofs, fail)
ifs.success.connect(npt.intake)
npt.success.connect(ofs.intake)
npt.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
