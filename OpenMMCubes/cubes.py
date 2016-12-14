import time
import numpy as np
from floe.api import OEMolComputeCube, parameter
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES
from OpenMMCubes.ports import OpenMMSystemOutput

from simtk import unit, openmm
from simtk.openmm import app

from openeye import oechem

import openmoltools
from openmoltools import forcefield_generators

# For parallel, import and inherit from ParallelOEMolComputeCube
class OpenMMComplexSetup(OEMolComputeCube):
    title = "Set up complex for simulation in OpenMM"
    description = """
    *Longform Description*
    Set up protein:ligand complex for simulation with OpenMM.
    """
    classification = [
        ["Testing", "OpenMM"],
        ["Testing", "Complex Setup"],
    ]
    tags = [tag for lists in classification for tag in lists]

    success = OpenMMSystemOutput("success")

    receptor = parameter.DataSetInputParameter(
        'receptor',
        required=True,
        help_text='Single Receptor to Dock Against')

    pH = parameter.DecimalParameter(
        'pH',
        default=7.4,
        help_text="Solvent pH used to select appropriate receptor protonation state.",
    )

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10,
        help_text="Padding around protein for solvent box (angstroms)",
    )

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=500,
        help_text="Salt concentration (millimolar)",
    )

    def begin(self):
        # Read receptor into molecule
        if in_orion():
            stream = StreamingDataset(self.args.receptor)
            mols = [ mol for mol in stream ]
            receptor = mols[0]
        else:
            receptor = oechem.OEMol()
            with oechem.oemolistream(self.args.receptor) as ifs:
                oechem.OEReadMolecule(ifs, receptor)

        # Write the receptor to a PDB
        pdbfilename = 'receptor.pdb'
        ofs = oechem.oemolostream(pdbfilename)
        oechem.OEWriteConstMolecule(ofs, receptor)
        ofs.close()

        # Read the PDB file into an OpenMM PDBFile object
        self.pdbfile = app.PDBFile(pdbfilename)

    def process(self, mol, port):
        try:
            # Create OpenMM Topology from OEMol
            molecule_topology = forcefield_generators.generateTopologyFromOEMol(mol)
            molecule_positions = unit.Quantity( np.zeros([mol.NumAtoms(),3], np.float64), unit.angstroms)
            for (index, atom) in enumerate(mol.GetAtoms()):
                [x,y,z] = mol.GetCoords(atom)
                molecule_positions[index,0] = x*unit.angstroms
                molecule_positions[index,1] = y*unit.angstroms
                molecule_positions[index,2] = z*unit.angstroms

            # Initialize a forcefield
            gaff_xml_filename = openmoltools.utils.get_data_filename("parameters/gaff.xml")
            forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml', gaff_xml_filename)
            forcefield.registerTemplateGenerator(forcefield_generators.gaffTemplateGenerator)

            # Set up protein
            modeller = app.Modeller(self.pdbfile.topology, self.pdbfile.positions)
            modeller.add(molecule_topology, molecule_positions)

            # Add missing hydrogens (if needed)
            modeller.addHydrogens(forcefield, self.args.pH)

            # Add solvent
            modeller.addSolvent(forcefield, model='tip3p', padding=self.args.solvent_padding,
                ionicStrength=unit.Quantity(self.args.salt_concentration, unit.millimolar)
                )

            # Generate System
            system = forcefield.createSystem(modeller.getTopology(), nonbondedMethod=app.PME, nonbondedCutoff=10.0*unit.angstroms, constraints=app.HBonds)

            # Emit the serialized system.
            self.success.emit(system)
        except Exception as e:
            # Attach error message to the molecule that failed
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)
