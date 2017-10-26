from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMminimizeCube, OpenMMnvtCube, OpenMMnptCube
from ComplexPrepCubes.cubes import HydrationCube, ComplexPrep, ForceFieldPrep
from ComplexPrepCubes.port import ProteinReader
from LigPrepCubes.ports import LigandReader
from LigPrepCubes.cubes import LigChargeCube
from YankCubes.cubes import  SyncBindingFECube, YankBindingFECube

job = WorkFloe('Yank Binding Affinity')

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

job.classification = [['BindingFreeEnergy', 'Yank']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = LigandReader("LigandReader", title="Ligand Reader")
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

# COMPLEX SETTING

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexPrep("Complex")

# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvateComplex = HydrationCube("HydrationComplex", title="HydrationComplex")

# Complex Force Field Application
ffComplex = ForceFieldPrep("ForceFieldComplex", title="ForceFieldComplex")
ffComplex.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='amber99sbildn.xml')
ffComplex.promote_parameter('solvent_forcefield', promoted_name='solvent_ff', default='tip3p.xml')
ffComplex.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='GAFF2')
ffComplex.promote_parameter('other_forcefield', promoted_name='other_ff', default='GAFF2')

# Minimization
minComplex = OpenMMminimizeCube('minComplex', title='MinimizeComplex')
minComplex.promote_parameter('restraints', promoted_name='m_restraints', default="noh (ligand or protein)",
                             description='Select mask to apply restarints')
minComplex.promote_parameter('restraintWt', promoted_name='m_restraintWt', default=5.0,
                             description='Restraint weight')

# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmupComplex = OpenMMnvtCube('warmupComplex', title='warmupComplex')
warmupComplex.promote_parameter('time', promoted_name='warm_psec', default=20.0,
                                description='Length of MD run in picoseconds')
warmupComplex.promote_parameter('restraints', promoted_name='w_restraints', default="noh (ligand or protein)",
                                description='Select mask to apply restarints')
warmupComplex.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.0, description='Restraint weight')
warmupComplex.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=0,
                                description='Trajectory saving interval')
warmupComplex.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=0,
                                description='Reporter saving interval')
warmupComplex.promote_parameter('outfname', promoted_name='w_outfname', default='warmup',
                                description='Equilibration suffix name')

# The system is equilibrated at the right pressure and temperature in 3 stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Equilibration stage 1
equil1Complex = OpenMMnptCube('equil1Complex', title='equil1Complex')
equil1Complex.promote_parameter('time', promoted_name='eq1_psec', default=20.0,
                                description='Length of MD run in picoseconds')
equil1Complex.promote_parameter('restraints', promoted_name='eq1_restraints', default="noh (ligand or protein)",
                                description='Select mask to apply restarints')
equil1Complex.promote_parameter('restraintWt', promoted_name='eq1_restraintWt', default=2.0, description='Restraint weight')
equil1Complex.promote_parameter('trajectory_interval', promoted_name='eq1_trajectory_interval', default=0,
                                description='Trajectory saving interval')
equil1Complex.promote_parameter('reporter_interval', promoted_name='eq1_reporter_interval', default=0,
                                description='Reporter saving interval')
equil1Complex.promote_parameter('outfname', promoted_name='eq1_outfname', default='equil1',
                                description='Equilibration suffix name')

# NPT Equilibration stage 2
equil2Complex = OpenMMnptCube('equil2Complex', title='equil2Complex')
equil2Complex.promote_parameter('time', promoted_name='eq2_psec', default=20.0,
                                description='Length of MD run in picoseconds')
equil2Complex.promote_parameter('restraints', promoted_name='eq2_restraints', default="noh (ligand or protein)",
                                description='Select mask to apply restarints')
equil2Complex.promote_parameter('restraintWt', promoted_name='eq2_restraintWt', default=0.5,
                                description='Restraint weight')
equil2Complex.promote_parameter('trajectory_interval', promoted_name='eq2_trajectory_interval', default=0,
                                description='Trajectory saving interval')
equil2Complex.promote_parameter('reporter_interval', promoted_name='eq2_reporter_interval', default=0,
                                description='Reporter saving interval')
equil2Complex.promote_parameter('outfname', promoted_name='eq2_outfname', default='equil2',
                                description='Equilibration suffix name')

# NPT Equilibration stage 3
equil3Complex = OpenMMnptCube('equil3Complex', title='equil3Complex')
equil3Complex.promote_parameter('time', promoted_name='eq3_psec', default=20.0,
                                description='Length of MD run in picoseconds')
equil3Complex.promote_parameter('restraints', promoted_name='eq3_restraints', default="ca_protein or (noh ligand)",
                                description='Select mask to apply restarints')
equil3Complex.promote_parameter('restraintWt', promoted_name='eq3_restraintWt', default=0.1,
                                description='Restraint weight')
equil3Complex.promote_parameter('trajectory_interval', promoted_name='eq3_trajectory_interval', default=0,
                                description='Trajectory saving interval')
equil3Complex.promote_parameter('reporter_interval', promoted_name='eq3_reporter_interval', default=0,
                                description='Reporter saving interval')
equil3Complex.promote_parameter('outfname', promoted_name='eq3_outfname', default='equil3',
                                description='Equilibration suffix name')

# LIGAND SETTING

# Solvate Ligands
solvateLigand = HydrationCube("HydrationLigand", title="HydrationLigand")
solvateLigand.promote_parameter('ref_structure', promoted_name='ref_structure', default='False')

# Ligand Force Field Application
ffLigand = ForceFieldPrep("ForceFieldLigand", title="ForceFieldLigand")
ffLigand.promote_parameter('solvent_forcefield', promoted_name='solvent_ff',
                           default=ffComplex.promoted_parameters['solvent_forcefield']['default'])
ffLigand.promote_parameter('ligand_forcefield', promoted_name='ligand_ff',
                           default=ffComplex.promoted_parameters['ligand_forcefield']['default'])
ffLigand.promote_parameter('other_forcefield', promoted_name='other_ff',
                           default=ffComplex.promoted_parameters['other_forcefield']['default'])

# Ligand Minimization
minimizeLigand = OpenMMminimizeCube("MinimizeLigand")
minimizeLigand.promote_parameter('restraints', promoted_name='m_restraints', default="noh ligand",
                                 description='Select mask to apply restraints')
minimizeLigand.promote_parameter('restraintWt', promoted_name='m_restraintWt', default=5.0,
                                 description='Restraint weight in kcal/(mol A^2')

# Ligand NVT Warm-up
warmupLigand = OpenMMnvtCube('warmupLigand', title='warmupLigand')
warmupLigand.promote_parameter('time', promoted_name='warm_psec', default=20.0,
                               description='Length of MD run in picoseconds')
warmupLigand.promote_parameter('restraints', promoted_name='w_restraints', default="noh ligand",
                               description='Select mask to apply restarints')
warmupLigand.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.0,
                               description='Restraint weight in kcal/(mol A^2')
warmupLigand.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=0,
                               description='Trajectory saving interval')
warmupLigand.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=0,
                               description='Reporter saving interval')
warmupLigand.promote_parameter('outfname', promoted_name='w_outfname', default='warmup',
                               description='Equilibration suffix name')
# warmupLigand.promote_parameter('center', promoted_name='center', default=True)

# Ligand NPT Equilibration stage
equilLigand = OpenMMnptCube('equilLigand', title='equilLigand')
equilLigand.promote_parameter('time', promoted_name='eq_psec', default=20.0,
                              description='Length of MD run in picoseconds')
equilLigand.promote_parameter('restraints', promoted_name='eq_restraints', default="noh ligand",
                              description='Select mask to apply restraints')
equilLigand.promote_parameter('restraintWt', promoted_name='eq_restraintWt', default=0.1,
                              description='Restraint weight in kcal/(mol A^2')
equilLigand.promote_parameter('trajectory_interval', promoted_name='eq_trajectory_interval', default=0,
                              description='Trajectory saving interval')
equilLigand.promote_parameter('reporter_interval', promoted_name='eq_reporter_interval', default=0,
                              description='Reporter saving interval')
equilLigand.promote_parameter('outfname', promoted_name='eq_outfname', default='equil',
                              description='Equilibration suffix name')

sync = SyncBindingFECube("SyncCube")

yank = YankBindingFECube("YankABFE")

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iprot, iligs, chargelig, complx, solvateComplex, ffComplex,
              minComplex, warmupComplex, equil1Complex, equil2Complex, equil3Complex,
              solvateLigand, ffLigand, minimizeLigand, warmupLigand, equilLigand,
              sync, yank, ofs, fail)

# Connections
iprot.success.connect(complx.system_port)
iligs.success.connect(chargelig.intake)
chargelig.success.connect(complx.intake)
# Complex Connections
complx.success.connect(solvateComplex.intake)
solvateComplex.success.connect(ffComplex.intake)
ffComplex.success.connect(minComplex.intake)
minComplex.success.connect(warmupComplex.intake)
warmupComplex.success.connect(equil1Complex.intake)
equil1Complex.success.connect(equil2Complex.intake)
equil2Complex.success.connect(equil3Complex.intake)
equil3Complex.success.connect(sync.intake)
# Ligand Connections
chargelig.success.connect(solvateLigand.intake)
solvateLigand.success.connect(ffLigand.intake)
ffLigand.success.connect(minimizeLigand.intake)
minimizeLigand.success.connect(warmupLigand.intake)
warmupLigand.success.connect(equilLigand.intake)
equilLigand.success.connect(sync.solvated_ligand_in_port)
# SYNC
sync.solvated_lig_complex_out_port.connect(yank.intake) 
yank.success.connect(ofs.intake)
yank.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
