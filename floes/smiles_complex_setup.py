from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import ChargeMCMol, SMIRFFParameterization, GAFFParameterization, FREDDocking

job = WorkFloe("SmilesComplexPrep")

job.description = """
Parse a SMILES string, dock
"""

job.classification = [
    ["Testing", "Complex Setup"],
]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="File containing SMILES")

charge = ChargeMCMol('charge')

fred = FREDDocking('fred')
fred.promote_parameter('receptor', promoted_name='receptor', description='Receptor OEB')

smirff = SMIRFFParameterization('smirff')
smirff.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="SMIRFF FFXML")

complex_setup = OpenMMComplexSetup("complex_setup")
complex_setup.promote_parameter('protein', promoted_name='protein', description="PDB of protein structure")
complex_setup.promote_parameter('pH', promoted_name='pH')
complex_setup.promote_parameter('solvent_padding', promoted_name='solvent_padding')
complex_setup.promote_parameter('salt_concentration', promoted_name='salt_conc')
complex_setup.promote_parameter('protein_forcefield', promoted_name='protein_ff')
complex_setup.promote_parameter('solvent_forcefield', promoted_name='solvent_ff')

ofs = OEMolOStreamCube('ofs')
ofs.set_parameters(data_out="complex.oeb.gz")

job.add_cubes(ifs, charge, fred, smirff, complex_setup, ofs)
ifs.success.connect(charge.intake)
charge.success.connect(fred.intake)
fred.success.connect(smirff.intake))
smirff.success.connect(complex_setup.intake)
complex_setup.success.connect(ofs.intake)
if __name__ == "__main__":
    job.run()
