from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from ComplexPrepCubes.cubes import Reader, Splitter, SolvationCube, LigChargeCube, ComplexPrep, ForceFieldPrep

job = WorkFloe("ComplexPrep")

job.description = """
Complex Preparation Workflow

Ex. python floes/openmm_complex_prep.py --protein protein.oeb
--Ligands-data_in ligands.oeb  --ofs-data_out complex.oeb

Parameters:
-----------
protein (file): OEB file of the prepared protein
ligands (file): OEB file of the prepared ligands


Outputs:
--------
ofs: Output file
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = OEMolIStreamCube("Ligands")
#iligs.promote_parameter("data_in", promoted_name="ligands", description="Input ligands")

chargelig = LigChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=5)

# Protein Setting
isys = Reader("ProteinReader")
isys.promote_parameter("data_in", promoted_name="protein", description="Protein file name")
isys.promote_parameter("protein_suffix", promoted_name="protein_suffix", default='MCL1', description="Protein suffix")

splitter = Splitter("Splitter")
solvate = SolvationCube("Solvation")

# Complex Setting
complx = ComplexPrep("Complex")
ff = ForceFieldPrep("ForceField")

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

job.add_cubes(isys, splitter, solvate, iligs, chargelig, complx, ff, ofs)

isys.success.connect(splitter.intake)
splitter.success.connect(solvate.intake)
solvate.success.connect(complx.system_port)
iligs.success.connect(chargelig.intake)
chargelig.success.connect(complx.intake)
complx.success.connect(ff.intake)
ff.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
