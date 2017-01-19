import time
import traceback
import numpy as np
from floe.api import OEMolComputeCube, parameter, MoleculeInputPort, BinaryMoleculeInputPort, BinaryOutputPort, OutputPort
from floe.api.orion import in_orion, StreamingDataset
from floe.constants import BYTES
from OpenMMCubes.ports import OpenMMSystemOutput, OpenMMSystemInput
from simtk import unit, openmm
from simtk.openmm import app

from openeye import oechem
import os, smarty, parmed, pdbfixer
from openmoltools import forcefield_generators
import lzma
from simtk.openmm import XmlSerializer

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

    #intake = MoleculeInputPort('intake')
    #success = OpenMMSystemOutput('success')
    #failure = OpenMMSystemOutput('failure')

    protein = parameter.DataSetInputParameter(
        'protein',
        required=True,
        help_text='Single protein to Dock Against')

    #ligand = parameter.DataSetInputParameter(
    #    'ligand',
    #    required=True,
    #    help_text='Docked ligands')

    pH = parameter.DecimalParameter(
        'pH',
        default=7.0,
        help_text="Solvent pH used to select appropriate protein protonation state.",
    )

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10,
        help_text="Padding around protein for solvent box (angstroms)",
    )

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=50,
        help_text="Salt concentration (millimolar)",
    )

    molecule_forcefield = parameter.DataSetInputParameter(
        'molecule_forcefield',
    #    required=True,
        help_text='Forcefield parameters for molecule'
    )

    protein_forcefield = parameter.DataSetInputParameter(
        'protein_forcefield',
    #    required=True,
        default='amber99sbildn.xml',
        help_text='Forcefield parameters for protein'
    )

    solvent_forcefield = parameter.DataSetInputParameter(
        'solvent_forcefield',
    #   required=True,
        default='tip3p.xml',
        help_text='Forcefield parameters for solvent'
    )

    def begin(self):
        pdbfilename = 'protein.pdb'

        # Write the protein to a PDB
        if in_orion():
            stream = StreamingDataset(self.args.protein, input_format=".pdb")
            stream.download_to_file(pdbfilename)
        else:
            protein = oechem.OEMol()
            with oechem.oemolistream(self.args.protein) as ifs:
                if not oechem.OEReadMolecule(ifs, protein):
                    raise RuntimeError("Error reading molecule")
                with oechem.oemolostream(pdbfilename) as ofs:
                    res = oechem.OEWriteMolecule(ofs, protein)
                    if res != oechem.OEWriteMolReturnCode_Success:
                        raise RuntimeError("Error writing protein: {}".format(res))

        # Read the PDB file into an OpenMM PDBFile object
        self.pdbfile = app.PDBFile(pdbfilename)

    def process(self, mol, port):
        try:
            # Generature ligand structure
            #ifs = oechem.oemolistream(self.args.ligand)
            #flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
            #ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
            #oechem.OEReadMolecule(ifs, mol)
            #oechem.OETriposAtomNames(mol)

            from smarty.forcefield import ForceField
            #mol_ff = ForceField(self.args.molecule_forcefield)
            molecule_forcefield = os.path.join('OpenMMCubes', 'tests', 'input', 'forcefield','smirff99Frosst.ffxml')
            mol_ff = ForceField(molecule_forcefield)
            mol_top, mol_sys, mol_pos = smarty.forcefield_utils.create_system_from_molecule(mol_ff, mol)
            molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys)

            #Alter molecule residue name for easy selection
            molecule_structure.residues[0].name = "MOL"

            #Generate protein Structure object
            forcefield = app.ForceField(self.args.protein_forcefield, self.args.solvent_forcefield)
            system = forcefield.createSystem( self.pdbfile.topology )
            protein_structure = parmed.openmm.load_topology( self.pdbfile.topology, system )

            # Merge structures
            structure = protein_structure + molecule_structure
            topology = structure.topology

            #Concatenate positions arrays
            positions_unit = unit.angstroms
            positions0_dimensionless = np.array( self.pdbfile.positions / positions_unit )
            positions1_dimensionless = np.array( mol_pos / positions_unit )

            coordinates = np.vstack((positions0_dimensionless,positions1_dimensionless))
            natoms = len(coordinates)
            positions = np.zeros([natoms,3], np.float32)
            for index in range(natoms):
                (x,y,z) = coordinates[index]
                positions[index,0] = x
                positions[index,1] = y
                positions[index,2] = z
            positions = unit.Quantity(positions, positions_unit)

            #Store in Structure object
            structure.coordinates = coordinates

            # Save to PDB
            structure.save('tmp.pdb', overwrite=True)

            # Solvate with PDBFixer
            fixer = pdbfixer.PDBFixer('tmp.pdb')
            #fixer.findMissingResidues()
            #fixer.findMissingAtoms()
            #fixer.addMissingAtoms()
            fixer.addMissingHydrogens(self.args.pH)
            fixer.addSolvent(padding=unit.Quantity(self.args.solvent_padding, unit.angstroms),
                            ionicStrength=unit.Quantity(self.args.salt_concentration, unit.millimolar)
                            )

            # Load PDBFixer object back to Structure
            struct = parmed.openmm.load_topology(fixer.topology, xyz=fixer.positions)

            # Remove ligand from protein Structure by AmberMask selection
            struct.strip(":MOL")

            # Regenerate openMM System to parameterize solvent
            system = forcefield.createSystem(struct.topology, rigidWater=False)

            # Regenerate parameterized protein structure
            protein_structure = parmed.openmm.load_topology( struct.topology, system=system,
                                                            xyz=struct.positions, box=struct.box )

            # Remerge with ligand structure
            combined_structure = protein_structure + molecule_structure

            # Restore initial positions and box dimensions
            combined_structure.positions = fixer.positions
            combined_structure.box = struct.box

            # Remove temporary pdb
            os.remove('tmp.pdb')
            combined_structure.save('combined_structure.pdb', overwrite=True)

            # Regenerate OpenMM system with parmed
            system = combined_structure.createSystem(nonbondedMethod=app.PME,
                                                    nonbondedCutoff=10.0*unit.angstroms,
                                                    constraints=app.HBonds)

            complex_mol = oechem.OEMol()
            output = OpenMMSystemOutput('output')
            with oechem.oemolistream('combined_structure.pdb') as ifs:
                if not oechem.OEReadMolecule(ifs, complex_mol):
                    raise RuntimeError("Error reading complex pdb")
                with oechem.oemolostream('setup.oeb.gz') as ofs:
                    complex_mol.SetData(oechem.OEGetTag('system'), output.encode(system))
                    self.success.emit(complex_mol)
                    #os.remove('combined_structure.pdb')

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

class OpenMMSimulation(OEMolComputeCube):
    title = "Run simulation in OpenMM"
    description = """
    *Longform Description*
    Run simulation with OpenMM for protein:ligand complex.
    """
    classification = [
        ["Testing", "OpenMM"],
        ["Testing", "Simulation"],
    ]
    tags = [tag for lists in classification for tag in lists]

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300,
        help_text="Temperature (Kelvin)"
    )

    steps = parameter.IntegerParameter(
        'steps',
        default=10000,
        help_text="Number of MD steps"
    )

    def begin(self):
        #pdbfilename = 'complex.pdb'        # Write the protein to a PDB
        #if in_orion():
        #    stream = StreamingDataset(self.args.complex_pdb, input_format=".pdb")
        #    stream.download_to_file(pdbfilename)
        #else:
        #    mol = oechem.OEMol()
        #    with oechem.oemolistream(self.args.complex_pdb) as ifs:
        #        if not oechem.OEReadMolecule(ifs, mol):
        #            raise RuntimeError("Error reading molecule")
        #        with oechem.oemolostream(pdbfilename) as ofs:
        #            res = oechem.OEWriteMolecule(ofs, mol)
        #            if res != oechem.OEWriteMolReturnCode_Success:
        #                raise RuntimeError("Error writing protein: {}".format(res))
        # Initialize openmm integrator
        self.integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)
        #self.pdbfile = app.PDBFile(pdbfilename)

    def process(self, mol, port):
        try:
            if not mol.HasData:
                raise RuntimeError('No Data found.')
            else:
                pdbfilename = 'complex.pdb'
                complex_mol = oechem.OEMol()
                with oechem.oemolostream(pdbfilename) as ofs:
                    oechem.OEWriteMolecule(ofs, mol)
                pdbfile = app.PDBFile(pdbfilename)
                serialized_system = mol.GetData(oechem.OEGetTag('system'))
                system = openmm.XmlSerializer.deserialize( serialized_system )

                # Initialize Simulation
                simulation = app.Simulation(pdbfile.getTopology(), system, self.integrator)

                try:
                    outfname = 'restart'
                    mol.GetData(oechem.OEGetTag('state'))
                    serialized_state = mol.GetData(oechem.OEGetTag('state'))
                    state = openmm.XmlSerializer.deserialize( serialized_state )
                    simulation.context.setState(state)
                except ValueError:
                    pass
                else:
                    outfname = 'simulation'
                    # Set initial positions and velocities then minimize
                    simulation.context.setPositions(pdbfile.getPositions())
                    simulation.context.setVelocitiesToTemperature(self.args.temperature*unit.kelvin)
                    simulation.minimizeEnergy()

            # Do MD simulation and report energies
            statereporter = app.StateDataReporter(outfname+'.log', 1000, step=True, potentialEnergy=True, temperature=True)
            simulation.reporters.append(statereporter)
            simulation.step(self.args.steps)

            # Save serialized State object
            state = simulation.context.getState( getPositions=True,
                                              getVelocities=True,
                                              getParameters=True )

            sim_mol = oechem.OEMol()
            output = OpenMMSystemOutput('output')
            with oechem.oemolistream(pdbfilename) as ifs:
                if not oechem.OEReadMolecule(ifs, sim_mol):
                    raise RuntimeError("Error reading {}".format(pdbfilename))
                with oechem.oemolostream(outfname+'.oeb.gz') as ofs:
                    sim_mol.SetData(oechem.OEGetTag('system'), output.encode(system))
                    sim_mol.SetData(oechem.OEGetTag('state'), output.encode(state))
                    with open(outfname+'.log', 'rb') as f:
                        sim_mol.SetData(oechem.OEGetTag('log'), f.read())
                    self.success.emit(sim_mol)

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(mol)
