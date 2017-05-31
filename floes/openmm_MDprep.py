from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import  OpenMMminimizeCube, OpenMMnvtCube, OpenMMnptCube
from ComplexPrepCubes.cubes import Reader, Splitter, SolvationCube, LigChargeCube, ComplexPrep, ForceFieldPrep


job = WorkFloe('Preparation MD')

job.description = """
Set up an OpenMM complex then minimize, warmup and equilibate MD

Ex. python floes/openmm_FrosstMD.py --ligand examples/data/toluene.pdb --protein examples/data/T4-protein.pdb

Parameters:
-----------
ligands (file): oeb file of ligand posed in the protein active site.
protein (file): oeb file of the protein structure, assumed to be pre-prepared

Optionals:
-----------
pH (float): Solvent pH used to select protein protonation states (default: 7.0)
solvent_padding (float): Padding around protein for solvent box (default: 10 angstroms)
salt_concentration (float): Salt concentration (default: 50 millimolar)
molecule_forcefield (file): Smarty parsable FFXML file containining parameters for the molecule (default: smirff99Frosst.ffxml)
protein_forcefield (file): XML file containing forcefield parameters for protein (default: amber99sbildn.xml)
solvent_forcefield (file): XML file containing forcefield parameter for solvent (default: tip3p.xml)

Outputs:
--------
ofs: Outputs a <IDTag>-complex.oeb.gz file containing the
OpenMM System, State, and ParmEd Structure of the protein:ligand complex,
packaged with the OEMol.
"""

job.classification = [['Complex Setup', 'FrosstMD']]
job.tags = [tag for lists in job.classification for tag in lists]


# Ligand setting
iligs = OEMolIStreamCube("Ligands")
# iligs.promote_parameter("data_in", promoted_name="ligand", description="PDB of docked ligand")
chargelig = LigChargeCube("LigCharge")

chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=5)

# Protein Setting
isys = Reader("ProteinReader")
isys.promote_parameter("data_in", promoted_name="protein", description="PDB Protein file name")
isys.promote_parameter("protein_suffix", promoted_name="protein_suffix", default='1HVL', description="Protein suffix")

# Splitter
splitter = Splitter("Splitter")

# Solvate
solvate = SolvationCube("Solvation")
solvate.promote_parameter('pH', promoted_name='pH')
solvate.promote_parameter('solvent_padding', promoted_name='solvent_padding')
solvate.promote_parameter('salt_concentration', promoted_name='salt_conc', default=100)

# Complex Setting
complx = ComplexPrep("Complex")

# Force Field Application
ff = ForceFieldPrep("ForceField")
ff.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='amber99sbildn.xml')
ff.promote_parameter('solvent_forcefield', promoted_name='solvent_ff', default='tip3p.xml')
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='GAFF2')

minComplex = OpenMMminimizeCube('minComplex', title='Minimize')
minComplex.promote_parameter('steps', promoted_name='steps', default=30000)
minComplex.promote_parameter('restraints', promoted_name='w_restraints', default="noh (ligand and protein)", description='Select mask to apply restarints')
minComplex.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=5.0, description='Restraint weight')

warmup = OpenMMnvtCube('warmup', title='warmup')
warmup.promote_parameter('time', promoted_name='warm_psec', default=2.0,  description='Length of MD run in picoseconds')
warmup.promote_parameter('restraints', promoted_name='w_restraints', default="noh (ligand and protein)", description='Select mask to apply restarints')
warmup.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.0, description='Restraint weight')
warmup.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=1000, description='Trajectory saving interval')
warmup.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=10000, description='Reporter saving interval')
warmup.promote_parameter('outfname', promoted_name='w_outfname', default='warmup', description='Equilibration suffix name')

equil1 = OpenMMnptCube('equil1', title='equil1')
equil1.promote_parameter('time', promoted_name='equil1_psec', default=2.0, description='Length of MD run in picoseconds')
equil1.promote_parameter('restraints', promoted_name='eq1_restraints', default="noh (ligand and protein)", description='Select mask to apply restarints')
equil1.promote_parameter('restraintWt', promoted_name='eq1_restraintWt', default=2.0, description='Restraint weight')
equil1.promote_parameter('trajectory_interval', promoted_name='eq_1trajectory_interval', default=1000, description='Trajectory saving interval')
equil1.promote_parameter('reporter_interval', promoted_name='eq1_reporter_interval', default=10000, description='Reporter saving interval')
equil1.promote_parameter('outfname', promoted_name='eq1_outfname', default='equil1', description='Equilibration suffix name')

equil2 = OpenMMnptCube('equil2', title='equil2')
equil2.promote_parameter('time', promoted_name='equil2_psec', default=2.0, description='Length of MD run in picoseconds')
equil2.promote_parameter('restraints', promoted_name='eq2_restraints', default="noh (ligand and protein)", description='Select mask to apply restarints')
equil2.promote_parameter('restraintWt', promoted_name='eq2_restraintWt', default=0.5, description='Restraint weight')
equil2.promote_parameter('trajectory_interval', promoted_name='eq2_trajectory_interval', default=1000, description='Trajectory saving interval')
equil2.promote_parameter('reporter_interval', promoted_name='eq2_reporter_interval', default=10000, description='Reporter saving interval')
equil2.promote_parameter('outfname', promoted_name='eq2_outfname', default='equil2', description='Equilibration suffix name')

equil3 = OpenMMnptCube('equil3', title='equil3')
equil3.promote_parameter('time', promoted_name='equil3_psec', default=2.0, description='Length of MD run in picoseconds')
equil3.promote_parameter('restraints', promoted_name='eq3_restraints', default="ca_protein and (noh ligand)", description='Select mask to apply restarints')
equil3.promote_parameter('restraintWt', promoted_name='eq3_restraintWt', default=0.1, description='Restraint weight')
equil3.promote_parameter('trajectory_interval', promoted_name='eq3_trajectory_interval', default=1000, description='Trajectory saving interval')
equil3.promote_parameter('reporter_interval', promoted_name='eq3_reporter_interval', default=10000, description='Reporter saving interval')
equil3.promote_parameter('outfname', promoted_name='eq3_outfname', default='equil3', description='Equilibration suffix name')


ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')


job.add_cubes(isys, splitter, solvate, iligs, chargelig, complx, ff,
              minComplex, warmup, equil1, equil2, equil3, ofs, fail)


isys.success.connect(splitter.intake)
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
