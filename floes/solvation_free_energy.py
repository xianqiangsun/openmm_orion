from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolOStreamCube
from ComplexPrepCubes.cubes import SolvationCube, ForceFieldPrep
from LigPrepCubes.cubes import LigChargeCube
from LigPrepCubes.ports import LigandReader
from YankCubes.cubes import YankSolvationFECube
from OpenMMCubes.cubes import OpenMMminimizeCube, OpenMMnvtCube, OpenMMnptCube

job = WorkFloe("SolvationFreeEnergy")

job.description = """
Solvation Free Energy Calculation of small molecules

Ex. python floes/solvation_free_energy --ligands ligands.oeb
--ofs-data_out fe.oeb

Parameters:
-----------
ligands (file): OEB file of the prepared ligands

Outputs:
--------
ofs: Output file
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = LigandReader("Ligands", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")

chargelig = LigChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)

solvate = SolvationCube("Solvation")
solvate.promote_parameter("density", promoted_name="density", title="Solution density in g/ml", default=1.0,
                          description="Solution Density in g/ml")
solvate.promote_parameter("solvents", promoted_name="solvents", title="Solvent components",
                          default='[H]O[H], ClC(Cl)Cl, CS(=O)C, c1ccccc1',
                          description="Comma separated smiles strings of solvent components")
solvate.promote_parameter("molar_fractions", promoted_name="molar_fractions",
                          title="Molar fractions",
                          default='1.0, 0.0, 0.0, 0.0',
                          description="Comma separated strings of solvent molar fractions")
solvate.promote_parameter('distance_between_atoms', promoted_name='distance_between_atoms', default=2.5)
solvate.promote_parameter("padding_distance", promoted_name="padding_distance", default=11.0,
                          description="The largest dimension (in A) of the solute (along the x, y, or z axis) "
                                      "is determined, and a cubic box of size "
                                      "(largest dimension)+2*padding is used")



ff = ForceFieldPrep("ForceField")
# ff.promote_parameter('ligand_forcefield', promoted_name='ligand_forcefield', default='SMIRNOFF')

# Minimization
minimize = OpenMMminimizeCube("Minimize")
minimize.promote_parameter('restraints', promoted_name='m_restraints', default="noh ligand",
                           description='Select mask to apply restarints')
minimize.promote_parameter('restraintWt', promoted_name='m_restraintWt', default=5.0,
                           description='Restraint weight in kcal/(mol A^2')

# NVT Warm-up
warmup = OpenMMnvtCube('warmup', title='warmup')
warmup.promote_parameter('time', promoted_name='warm_psec', default=20.0,
                         description='Length of MD run in picoseconds')
warmup.promote_parameter('restraints', promoted_name='w_restraints', default="noh ligand",
                         description='Select mask to apply restarints')
warmup.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.0,
                         description='Restraint weight in kcal/(mol A^2')
warmup.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=0,
                         description='Trajectory saving interval')
warmup.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=0,
                         description='Reporter saving interval')
warmup.promote_parameter('outfname', promoted_name='w_outfname', default='warmup',
                         description='Equilibration suffix name')
warmup.promote_parameter('center', promoted_name='center', default=True)

# NPT Equilibration stage
equil = OpenMMnptCube('equil', title='equil')
equil.promote_parameter('time', promoted_name='eq_psec', default=20.0,
                        description='Length of MD run in picoseconds')
equil.promote_parameter('restraints', promoted_name='eq_restraints', default="noh ligand",
                        description='Select mask to apply restarints')
equil.promote_parameter('restraintWt', promoted_name='eq_restraintWt', default=0.1,
                        description='Restraint weight in kcal/(mol A^2')
equil.promote_parameter('trajectory_interval', promoted_name='eq_trajectory_interval', default=0,
                        description='Trajectory saving interval')
equil.promote_parameter('reporter_interval', promoted_name='eq_reporter_interval', default=0,
                        description='Reporter saving interval')
equil.promote_parameter('outfname', promoted_name='eq_outfname', default='equil',
                        description='Equilibration suffix name')

solvationfe = YankSolvationFECube("SovationFE")
solvationfe.promote_parameter('iterations', promoted_name='iterations', default=1000)
solvationfe.promote_parameter('nonbondedCutoff', promoted_name='nonbondedCutoff', default=10.0)

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iligs, chargelig, solvate, ff, minimize, warmup, equil, solvationfe, ofs, fail)

iligs.success.connect(chargelig.intake)
chargelig.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(minimize.intake)
minimize.success.connect(warmup.intake)
warmup.success.connect(equil.intake)
equil.success.connect(solvationfe.intake)
solvationfe.success.connect(ofs.intake)
solvationfe.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()