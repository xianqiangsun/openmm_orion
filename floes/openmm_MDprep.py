from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMminimizeCube, OpenMMnvtCube, OpenMMnptCube
from ComplexPrepCubes.cubes import ProteinReader, Splitter, SolvationCube, \
    ComplexPrep, ForceFieldPrep

from LigPrepCubes.cubes import LigChargeCube

job = WorkFloe('Preparation MD')

job.description = """
Set up an OpenMM complex then minimize, warm up and equilibrate a system by using three equilibration stages

Ex: python floes/openmm_MDprep.py --ligands ligands.oeb --protein protein.oeb --ofs-data_out prep.oeb

Parameters:
-----------
ligands (file): oeb file of ligand posed in the protein active site.
protein (file): oeb file of the protein structure, assumed to be pre-prepared

Optionals:
-----------

Outputs:
--------
ofs: Outputs a ready system to MD production run
"""

job.classification = [['Complex Setup', 'FrosstMD']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand reading cube setting
iligs = OEMolIStreamCube("Ligands", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")


chargelig = LigChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)

# Protein Reading cube. The protein prefix parameter is used to select a name for the
# output system files
iprot = ProteinReader("ProteinReader")
iprot.promote_parameter("data_in", promoted_name="protein", title='Protein Input File',
                        description="Protein file name")
iprot.promote_parameter("protein_prefix", promoted_name="protein_prefix",
                        default='PRT', description="Protein prefix")

# The spitter cube will be used to pre-process the read in protein system
splitter = Splitter("Splitter")

# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvate = SolvationCube("Solvation")
solvate.promote_parameter('solvent_padding', promoted_name='solvent_padding', default=10)
solvate.promote_parameter('salt_concentration', promoted_name='salt_conc', default=100)

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexPrep("Complex")

# Force Field Application
ff = ForceFieldPrep("ForceField")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='amber99sbildn.xml')
ff.promote_parameter('solvent_forcefield', promoted_name='solvent_ff', default='tip3p.xml')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='GAFF2')
ff.promote_parameter('other_forcefield', promoted_name='other_ff', default='GAFF2')

# Minimization
minComplex = OpenMMminimizeCube('minComplex', title='Minimize')
minComplex.promote_parameter('steps', promoted_name='min_steps', default=30000)
minComplex.promote_parameter('restraints', promoted_name='m_restraints', default="noh (ligand or protein)",
                             description='Select mask to apply restarints')
minComplex.promote_parameter('restraintWt', promoted_name='m_restraintWt', default=5.0,
                             description='Restraint weight')

# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmup = OpenMMnvtCube('warmup', title='warmup')
warmup.promote_parameter('time', promoted_name='warm_psec', default=100.0,
                         description='Length of MD run in picoseconds')
warmup.promote_parameter('restraints', promoted_name='w_restraints', default="noh (ligand or protein)",
                         description='Select mask to apply restarints')
warmup.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.0, description='Restraint weight')
warmup.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=1000,
                         description='Trajectory saving interval')
warmup.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=10000,
                         description='Reporter saving interval')
warmup.promote_parameter('outfname', promoted_name='w_outfname', default='warmup',
                         description='Equilibration suffix name')

# The system is equilibrated at the right pressure and temperature in 3 stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Equilibration stage 1
equil1 = OpenMMnptCube('equil1', title='equil1')
equil1.promote_parameter('time', promoted_name='eq1_psec', default=100.0,
                         description='Length of MD run in picoseconds')
equil1.promote_parameter('restraints', promoted_name='eq1_restraints', default="noh (ligand or protein)",
                         description='Select mask to apply restarints')
equil1.promote_parameter('restraintWt', promoted_name='eq1_restraintWt', default=2.0, description='Restraint weight')
equil1.promote_parameter('trajectory_interval', promoted_name='eq1_trajectory_interval', default=1000,
                         description='Trajectory saving interval')
equil1.promote_parameter('reporter_interval', promoted_name='eq1_reporter_interval', default=10000,
                         description='Reporter saving interval')
equil1.promote_parameter('outfname', promoted_name='eq1_outfname', default='equil1',
                         description='Equilibration suffix name')

# NPT Equilibration stage 2
equil2 = OpenMMnptCube('equil2', title='equil2')
equil2.promote_parameter('time', promoted_name='eq2_psec', default=100.0,
                         description='Length of MD run in picoseconds')
equil2.promote_parameter('restraints', promoted_name='eq2_restraints', default="noh (ligand or protein)",
                         description='Select mask to apply restarints')
equil2.promote_parameter('restraintWt', promoted_name='eq2_restraintWt', default=0.5,
                         description='Restraint weight')
equil2.promote_parameter('trajectory_interval', promoted_name='eq2_trajectory_interval', default=1000,
                         description='Trajectory saving interval')
equil2.promote_parameter('reporter_interval', promoted_name='eq2_reporter_interval', default=10000,
                         description='Reporter saving interval')
equil2.promote_parameter('outfname', promoted_name='eq2_outfname', default='equil2',
                         description='Equilibration suffix name')

# NPT Equilibration stage 3
equil3 = OpenMMnptCube('equil3', title='equil3')
equil3.promote_parameter('time', promoted_name='eq3_psec', default=200.0,
                         description='Length of MD run in picoseconds')
equil3.promote_parameter('restraints', promoted_name='eq3_restraints', default="ca_protein or (noh ligand)",
                         description='Select mask to apply restarints')
equil3.promote_parameter('restraintWt', promoted_name='eq3_restraintWt', default=0.1,
                         description='Restraint weight')
equil3.promote_parameter('trajectory_interval', promoted_name='eq3_trajectory_interval', default=1000,
                         description='Trajectory saving interval')
equil3.promote_parameter('reporter_interval', promoted_name='eq3_reporter_interval', default=10000,
                         description='Reporter saving interval')
equil3.promote_parameter('outfname', promoted_name='eq3_outfname', default='equil3',
                         description='Equilibration suffix name')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iprot, splitter, solvate, iligs, chargelig, complx, ff,
              minComplex, warmup, equil1, equil2, equil3, ofs, fail)

iprot.success.connect(splitter.intake)
splitter.success.connect(solvate.intake)
solvate.success.connect(complx.system_port)
iligs.success.connect(chargelig.intake)
chargelig.success.connect(complx.intake)
complx.success.connect(ff.intake)
ff.success.connect(minComplex.intake)
minComplex.success.connect(warmup.intake)
warmup.success.connect(equil1.intake)
equil1.success.connect(equil2.intake)
equil2.success.connect(equil3.intake)
equil3.success.connect(ofs.intake)
equil3.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
