from __future__ import unicode_literals
from floe.api import WorkFloe, OEMolIStreamCube, OEMolOStreamCube, FileOutputCube, DataSetInputParameter, FileInputCube
from OpenMMCubes.cubes import OpenMMComplexSetup, OpenMMSimulation
from LigPrepCubes.cubes import ChargeMCMol, SMIRNOFFParameterization, GAFFParameterization, FREDDocking

job = WorkFloe("SmilesLigPrep")

job.description = """
Parse SMILES, generate multiconformer OEMolecules, assign partial charges, dock using
the FRED docking engine and parameterize the ligand with the SMIRNOFF forcefield parameters
"""

job.classification = [["Ligand Preparation"]]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = OEMolIStreamCube("ifs")
ifs.promote_parameter("data_in", promoted_name="ligand", description="File containing SMILES")

charge = ChargeMCMol('charge')

fred = FREDDocking('fred')
fred.promote_parameter('receptor', promoted_name='receptor', description='Receptor OEB')

smirnoff = smirnoffParameterization('smirnoff')
smirnoff.promote_parameter('molecule_forcefield', promoted_name='ffxml', description="smirnoff FFXML")

ofs = OEMolOStreamCube('ofs')
ofs.set_parameters(data_out="smirnoff.oeb.gz")

job.add_cubes(ifs, charge, fred, smirnoff, ofs)
ifs.success.connect(charge.intake)
charge.success.connect(fred.intake)
fred.success.connect(smirnoff.intake)
smirnoff.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
