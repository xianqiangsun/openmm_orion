from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMminimizeCube

job = WorkFloe("MDminimize")

job.description = """
Minimize an OpenMM-ready solvated complex

Ex: python floes/openmm_prepMDminimize.py --complex complex.oeb --ofs-data_out min.oeb --steps 1000`

Parameters:
-----------
complex (file): OEB file of the prepared protein:ligand complex

Optional:
--------
steps (int): Number of MD steps to minimize the system. If 0 until convergence will be reached

Outputs:
--------
ofs: Outputs the minimized system
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("complex", title="Complex Reader")
ifs.promote_parameter("data_in", promoted_name="complex", title='Complex Input File',
                      description="protein:ligand complex input file")

minComplex = OpenMMminimizeCube('minComplex')
minComplex.promote_parameter('steps', promoted_name='steps')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(ifs, minComplex, ofs, fail)
ifs.success.connect(minComplex.intake)
minComplex.success.connect(ofs.intake)
minComplex.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
