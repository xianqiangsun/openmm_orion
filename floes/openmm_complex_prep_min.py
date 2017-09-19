from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolOStreamCube
from ComplexPrepCubes.cubes import SolvationCube, ComplexPrep, ForceFieldPrep
from ComplexPrepCubes.port import ProteinReader
from OpenMMCubes.cubes import OpenMMminimizeCube
from LigPrepCubes.cubes import LigChargeCube
from LigPrepCubes.ports import LigandReader

job = WorkFloe("ComplexPrepMin")

job.description = """
Complex Preparation and Minimization Workflow

Ex. python floes/openmm_complex_prep.py --protein protein.oeb
--ligands ligands.oeb  --ofs-data_out min.oeb

Parameters:
-----------
protein (file): OEB file of the prepared protein
ligands (file): OEB file of the prepared ligands


Outputs:
--------
ofs: Output file of the assembled and minimized complex
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = LigandReader("Ligands", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")

chargelig = LigChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)

# Protein Setting
iprot = ProteinReader("ProteinReader")
iprot.promote_parameter("data_in", promoted_name="protein", title="Protein Input File", description="Protein file name")
iprot.promote_parameter("protein_prefix", promoted_name="protein_prefix", default='PRT',
                        description="Protein Prefix")

solvate = SolvationCube("Solvation")

# Complex Setting
complx = ComplexPrep("Complex")
ff = ForceFieldPrep("ForceField")

# Minimization
minComplex = OpenMMminimizeCube('minComplex')
minComplex.promote_parameter('steps', promoted_name='steps')

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iprot, solvate, iligs, chargelig, complx, ff, minComplex, ofs, fail)

iprot.success.connect(solvate.intake)
solvate.success.connect(complx.system_port)
iligs.success.connect(chargelig.intake)
chargelig.success.connect(complx.intake)
complx.success.connect(ff.intake)
ff.success.connect(minComplex.intake)
minComplex.success.connect(ofs.intake)
minComplex.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
