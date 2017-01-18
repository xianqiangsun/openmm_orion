from __future__ import unicode_literals
"""
Copyright (C) 2016 OpenEye Scientific Software
"""
from floe.api import WorkFloe, OEMolIStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from OpenMMCubes.ports import OpenMMSystemInput

job = WorkFloe("simulation_resume")

job.description = """
**Set up OpenMM complex for simulation**

Check out the awesome stuff at the [OpenMM website](http://openmm.org)
"""

job.classification = [
    ["OpenEye", "OpenMM"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="complex_pdb", description="complex pdb file")

md_sim = OpenMMSimulation('md_sim')
md_sim.promote_parameter('state', promoted_name='state')
md_sim.promote_parameter('system', promoted_name='system')
md_sim.promote_parameter('complex_pdb', promoted_name='complex_pdb')

#state_save = FileOutputCube('state_save')
#state_save.set_parameters(name="state_cont.xml.xz")

ofs = FileOutputCube('ofs')
ofs.set_parameters(name='chkpt_out.log')

job.add_cubes(ifs, md_sim, ofs)
ifs.success.connect(md_sim.intake)
md_sim.success.connect(ofs.intake)







if __name__ == "__main__":
    job.run()
