from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from OpenMMCubes.cubes import OpenMMComplexSetup
from LigPrepCubes.cubes import ChargeMCMol, LigandParameterization
from OpenMMCubes.cubes import  OpenMMminimizeCube, OpenMMnvtCube, OpenMMnptCube

job = WorkFloe('Preparation MD')

job.description = """
**Set up an OpenMM complex then minimize, warmup and equilibate MD **

This floe will do the following in each cube:
  (1) ifs: Read in the ligand file,
  (2) readin the Protein
  (3) Gaff: Parameterize the molecule with gaff and 
            generate the ParmEd Structure and attach it to the OEMol.
  (4) complex_setup: Paramterize the protein and merge with the molecule Structure,
        Using PDBFixer: add missing atoms, add hydrogens given a pH, and solvate with TIP3P.
        Attach tagged data containing the <IDTag>, <Structure> and <System>
  (5) Run minimization, warmup, and equilibration
  (5) ofs: Write out the OEMOl of the complex to a <IDTag>-complex.oeb.gz

Ex. `python floes/openmm_FrosstMD.py --ligand examples/data/toluene.pdb --protein examples/data/T4-protein.pdb`

Parameters:
-----------
ligand (file): PDB file of ligand posed in the protein active site.
protein (file): PDB file of the protein structure, *assumed to be `pre-prepared`*

*Optionals:
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

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="PDB of docked ligand")

charge = ChargeMCMol('charge')
charge.promote_parameter('max_conformers', promoted_name='max_conformers', description="Set the max number of conformers per ligand", default=1)
charge.promote_parameter('keep_conformers', promoted_name='keep_conformers', description="Set the number of conformers to keep")

lig_param = LigandParameterization('lig_param')
lig_param.promote_parameter('molecule_forcefield', promoted_name='molecule_forcefield', description='Forcefield for molecule')

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein', description="PDB of protein structure")
complex_setup.promote_parameter('pH', promoted_name='pH')
complex_setup.promote_parameter('solvent_padding', promoted_name='solvent_padding')
complex_setup.promote_parameter('salt_concentration', promoted_name='salt_conc', default=100)
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_ff')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_ff')


minComplex = OpenMMminimizeCube('minComplex', title='Minimize')
minComplex.promote_parameter('steps', promoted_name='steps')

warmup = OpenMMnvtCube('warmup', title='warmup')
warmup.promote_parameter('time', promoted_name='warm_psec', default=1.0)

equil = OpenMMnptCube('equil', title='equil')
equil.promote_parameter('time', promoted_name='equil1_psec', default=1.0,
      description='Length of MD run in picoseconds')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(ifs, charge, lig_param, complex_setup,
              minComplex, warmup, equil, ofs, fail)

ifs.success.connect(charge.intake)
charge.success.connect(lig_param.intake)
lig_param.success.connect(complex_setup.intake)
complex_setup.success.connect(minComplex.intake)
minComplex.success.connect(warmup.intake)
warmup.success.connect(equil.intake)
equil.success.connect(ofs.intake)
equil.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
