from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMnptCube

job = WorkFloe("NPT Run")

job.description = """
NPT simulation of an OpenMM-ready System

Ex: python floes/openmm_MDnpt.py --system complex.oeb --picosec 10

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
picosec (float): Number of picoseconds to warm up the complex
temperature (decimal): target final temperature in K
pressure (decimal): target final pressure in atm

Outputs:
--------
ofs: Outputs the constant temperature and pressure system
"""

job.classification = [['NPT']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

npt = OpenMMnptCube('npt')
npt.promote_parameter('time', promoted_name='picosec', default=10.0,
                      description='Length of MD run in picoseconds')
npt.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                      description='Selected temperature in K')
npt.promote_parameter('pressure', promoted_name='pressure', default=1.0,
                      description='Selected pressure in atm')

# Restraints
npt.promote_parameter('restraints', promoted_name='restraints', default="ca_protein or (noh ligand)",
                      description='Select mask to apply restraints')
npt.promote_parameter('restraintWt', promoted_name='restraintWt', default=5.0, description='Restraint weight')

# Trajectory and logging info frequency intervals
npt.promote_parameter('trajectory_interval', promoted_name='trajectory_interval', default=100,
                      description='Trajectory saving interval')
npt.promote_parameter('reporter_interval', promoted_name='reporter_interval', default=1000,
                      description='Reporter saving interval')
npt.promote_parameter('outfname', promoted_name='suffix', default='npt',
                      description='Equilibration suffix name')

npt.promote_parameter('tar', promoted_name='tar', default=True)


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
