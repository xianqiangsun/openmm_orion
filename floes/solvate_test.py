from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolOStreamCube
from ComplexPrepCubes.cubes import SolvationCube
from LigPrepCubes.ports import LigandReader

job = WorkFloe("Testing Solvation Cube")

ifs = LigandReader("solute", title="Solute Reader")
ifs.promote_parameter("data_in", promoted_name="solute", title='Solute Input File',
                      description="Solute input file")

ifs.promote_parameter("IDTag", promoted_name="IDTag", default=True)

solv = SolvationCube("SolvationCube")
solv.promote_parameter('solvents', promoted_name='solvents', default='[H]O[H], ClC(Cl)Cl, CS(=O)C, c1ccccc1, CCO, [Cl-]')
solv.promote_parameter('molar_fractions', promoted_name='molar_fractions', default='1.0,0.0,0.0,0.0,0.0,0.0')


#solv.promote_parameter("geometry", promoted_name="geometry", default="sphere")


ofs = OEMolOStreamCube('ofs', title='OFS-Success')
ofs.set_parameters(backend='s3')
fail = OEMolOStreamCube('fail', title='OFS-Failure')
fail.set_parameters(backend='s3')
fail.set_parameters(data_out='fail.oeb.gz')


job.add_cubes(ifs, solv, ofs, fail)

ifs.success.connect(solv.intake)
solv.success.connect(ofs.intake)
solv.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
