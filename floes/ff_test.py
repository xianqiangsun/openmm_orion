from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolOStreamCube
from LigPrepCubes.cubes import LigChargeCube
from ComplexPrepCubes.cubes import LigandReader, ProteinReader
from ComplexPrepCubes.cubes import ForceFieldPrep, SolvationCube

job = WorkFloe("Testing Force Field Application")

# ifs = LigandReader("ligand", title="Ligand Reader")
# ifs.promote_parameter("data_in", promoted_name="ligands", title='Ligand Input File',
#                       description="Ligand input file")
# ifs.promote_parameter("ID", promoted_name="ID", default=True)


ifs = ProteinReader("Protein Reader")
ifs.promote_parameter("data_in", promoted_name="protein", title='Protein Input File',
                      description="Protein input file")

ligchg = LigChargeCube("Ligand Charge")

solv = SolvationCube("Solvation")

ff = ForceFieldPrep("ForceField")

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

# job.add_cubes(ifs, ligchg, solv, ff, ofs, fail)
# ifs.success.connect(ligchg.intake)
# ligchg.success.connect(solv.intake)
# solv.success.connect(ff.intake)
# ff.success.connect(ofs.intake)
# ff.failure.connect(fail.intake)


job.add_cubes(ifs, solv, ff, ofs, fail)

ifs.success.connect(solv.intake)
solv.success.connect(ff.intake)
ff.success.connect(ofs.intake)
ff.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
