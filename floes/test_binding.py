from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube
from YankCubes.cubes import YankBindingFECube, SyncBindingFECube

job = WorkFloe("YANK small molecule absolute binding free energies")

job.description = """
# Compute small molecule absolute binding free energies using YANK.

This Floe processes the provided molecules, computes free energies of binding to a 
specified receptor, and appends the following SDData properties to the original molecules:
"""

iligand = OEMolIStreamCube("ReadingLigand")
iligand.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File",
                          description="Ligand file name")

icomplex = OEMolIStreamCube("ReadingComplex")
icomplex.promote_parameter("data_in", promoted_name="complex", title="Complex Input File",
                           description="complex file name")

sync = SyncBindingFECube("SyncCube")
yankabfe = YankBindingFECube("YankBindingFE")

ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')

fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iligand, icomplex, sync, yankabfe, ofs, fail)

icomplex.success.connect(sync.intake)
iligand.success.connect(sync.solvated_ligand_in_port)
sync.solvated_lig_complex_out_port.connect(yankabfe.intake)

yankabfe.success.connect(ofs.intake)
yankabfe.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()