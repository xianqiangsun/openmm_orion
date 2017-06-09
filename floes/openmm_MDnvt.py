from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMnvtCube

job = WorkFloe("NVT")

job.description = """
NVT simulation of an OpenMM-ready solvated complex

Ex: python floes/openmm_MDnvt.py --complex complex.oeb --ofs-data_out nvt.oeb --picosec 10.0

Parameters:
-----------
complex (file): OEB file of the prepared protein:ligand complex

Optional:
--------
picosec (float): Number of picoseconds to warm up the complex
temperature (decimal): target final temperature after warming

Outputs:
--------
ofs: Outputs the constant temperature and volume system
"""

job.classification = [['NVT']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")

nvt = OpenMMnvtCube('nvt')
nvt.promote_parameter('time', promoted_name='picosec', default=10.0)
nvt.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                      description='Selected temperature in K')
# Restraints
nvt.promote_parameter('restraints', promoted_name='restraints', default='noh (ligand or protein)')
nvt.promote_parameter('restraintWt', promoted_name='restraintWt', default=2.0)
# Trajectory and logging info frequency intervals
nvt.promote_parameter('trajectory_interval', promoted_name='trajectory_interval', default=10)
nvt.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=100,
                      description='Reporter saving interval')

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
