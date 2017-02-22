from __future__ import unicode_literals
"""
Written by John D. Chodera
"""
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from YankCubes.cubes import YankHydrationCube

job = WorkFloe("Yank hydration free energy")

job.description = """
Compute hydration free energies in explicit solvent using YANK.

Ex. `python floes/yank_hydration.py --molecules examples/data/freesolv_mini.oeb.gz --output output.sdf`

Parameters:
-----------
molecules (file): input molecules

*Optionals:
-----------
temperature (float): Temperature in Kelvin (default: 300.0)
pressure (float): Pressure in atmospheres (default: 1.0)
simulation_time (float): Simulation time in ns/replica (default: 0.010)
timestep (float): Timestep in fs (default: 2.0)
yaml_template (str): YANK YAML file to use (default template provided)

Outputs:
--------
Attaches to OEMols 'DeltaG_yank_hydration' and 'dDeltaG_yank_hydration'
"""

job.classification =[["YANK", "Hydration"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="molecules", description="Input molecules")

yank_cube = YankHydrationCube('yank_hydration')
yank_cube.promote_parameter('temperature', promoted_name='temperature')
yank_cube.promote_parameter('pressure', promoted_name='pressure')
yank_cube.promote_parameter('simulation_time', promoted_name='simulation_time')
yank_cube.promote_parameter('timestep', promoted_name='timestep')
yank_cube.promote_parameter('yaml_template', promoted_name='yaml_template')

ofs = OEMolOStreamCube("ofs")
ofs.promote_parameter("data_out", promoted_name="output", description="Output molecules")

cubes = [ifs, yank_cube, ofs]

job.add_cubes(*cubes)

ifs.success.connect(yank_cube.intake)
yank_cube.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
