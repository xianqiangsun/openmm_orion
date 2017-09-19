from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolOStreamCube
from ComplexPrepCubes.cubes import SolvationCube, ForceFieldPrep
from LigPrepCubes.cubes import LigChargeCube
from LigPrepCubes.ports import LigandReader
from YankCubes.cubes import YankSolvationCube

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
solvate.promote_parameter("solvents", promoted_name="solvents", title="Solvent components", default='[H]O[H], ClC(Cl)Cl, CS(=O)C, c1ccccc1',
                          description="Comma separated smiles strings of solvent components")
solvate.promote_parameter("molar_fractions", promoted_name="molar_fractions",
                          title="Molar fractions", default='1.0, 0.0, 0.0, 0.0',
                          description="Comma separated  strings of solvent molar fractions")

ff = ForceFieldPrep("ForceField")

solvationfe = YankSolvationCube("SovationFE")

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iligs, chargelig, solvate, ff,  solvationfe, ofs, fail)

iligs.success.connect(chargelig.intake)
chargelig.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(solvationfe.intake)
solvationfe.success.connect(ofs.intake)
solvationfe.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()