from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import ChargeMCMol, LigandParameterization, FREDDocking

job = WorkFloe("RunOpenMMSimulation")

job.description = """
**Run an OpenMM Simulation (default 500K steps / 1ns)**
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex", description="OEB of the protein:ligand complex")

md = OpenMMSimulation('md')
md.promote_parameter('steps', promoted_name='steps')
md.set_parameters(item_count=1)
md.set_parameters(prefetch_count=1)
md.set_parameters(item_timeout=28800) #8hrs / 1ns


ofs = OEMolOStreamCube('ofs')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail')
fail.set_parameters(data_out="sim_failed.oeb.gz")
fail.set_parameters(backend='s3')

job.add_cubes(ifs, md, ofs, fail)
ifs.success.connect(md.intake)
md.success.connect(ofs.intake)
md.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
