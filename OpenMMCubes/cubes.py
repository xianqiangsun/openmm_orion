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

    protein = parameter.DataSetInputParameter(
        'protein',
        required=True,
        help_text='Single protein to Dock Against')

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
        self.proteinpdb = app.PDBFile(pdbfilename)

    def process(self, mol, port):
        try:
            # Generate smarty ligand structure
            from smarty.forcefield import ForceField
            ffxml = mol.GetData(oechem.OEGetTag('forcefield')).encode()
            with open('mol_parameters.ffxml', 'wb') as out:
                out.write(ffxml)
            mol_ff = ForceField(open('mol_parameters.ffxml'))
            mol_top, mol_sys, mol_pos = smarty.forcefield_utils.create_system_from_molecule(mol_ff, mol)
            molecule_structure = parmed.openmm.load_topology(mol_top, mol_sys, xyz=mol_pos)
            #Alter molecule residue name for easy selection
            molecule_structure.residues[0].name = "MOL"

            #Generate protein Structure object
            forcefield = app.ForceField(self.args.protein_forcefield, self.args.solvent_forcefield)
            protein_system = forcefield.createSystem( self.proteinpdb.topology )
            protein_structure = parmed.openmm.load_topology( self.proteinpdb.topology,
                                                             protein_system,
                                                             xyz=self.proteinpdb.positions )

            # Merge structures
            pl_structure = protein_structure + molecule_structure

            #Concatenate positions arrays (ensures same units)
            positions_unit = unit.angstroms
            positions0_dimensionless = np.array( self.proteinpdb.positions / positions_unit )
            positions1_dimensionless = np.array( molecule_structure.positions / positions_unit )
            coordinates = np.vstack((positions0_dimensionless,positions1_dimensionless))
            natoms = len(coordinates)
            positions = np.zeros([natoms,3], np.float32)
            for index in range(natoms):
                (x,y,z) = coordinates[index]
                positions[index,0] = x
                positions[index,1] = y
                positions[index,2] = z
            positions = unit.Quantity(positions, positions_unit)

            #Store Structure object
            #structure.coordinates = coordinates
            pl_structure.positions = positions
            pl_structure.save('pl_tmp.pdb', overwrite=True)

            # Solvate with PDBFixer
            fixer = pdbfixer.PDBFixer('pl_tmp.pdb')
            #fixer.findMissingResidues()
            #fixer.findMissingAtoms()
            #fixer.addMissingAtoms()
            fixer.addMissingHydrogens(self.args.pH)
            fixer.addSolvent(padding=unit.Quantity(self.args.solvent_padding, unit.angstroms),
                            ionicStrength=unit.Quantity(self.args.salt_concentration, unit.millimolar)
                            )

            # Load PDBFixer object back to Structure
            tmp = parmed.openmm.load_topology(fixer.topology, xyz=fixer.positions)
            #Store positions, topology, and box vectors for solvated system
            full_positions = tmp.positions
            full_topology = tmp.topology
            full_box = tmp.box

            # Remove ligand from protein Structure by AmberMask selection
            tmp.strip(":MOL")
            tmp.save('nomol_tmp.pdb', overwrite=True)
            nomol = parmed.load_file('nomol_tmp.pdb')

            # Regenerate openMM System to parameterize solvent
            nomol_system = forcefield.createSystem(nomol.topology, rigidWater=False)

            # Regenerate parameterized solvated protein structure
            solv_structure = parmed.openmm.load_topology( nomol.topology,
                                                            nomol_system,
                                                            xyz=nomol.positions,
                                                            box=nomol.box )

            # Remerge with ligand structure
            full_structure = solv_structure + molecule_structure
            # Restore box dimensions
            full_structure.box = nomol.box
            # Save full structure
            full_structure.save('complex.pdb', overwrite=True)

            # Regenerate OpenMM system with parmed
            system = full_structure.createSystem(nonbondedMethod=app.PME,
                                                nonbondedCutoff=10.0*unit.angstroms,
                                                constraints=app.HBonds)

            self.log.info('Generated OpenMM system')
            self.log.info('Saving System to complex.oeb.gz')
            # Pack mol into oeb and emit system
            complex_mol = oechem.OEMol()
            output = OpenMMSystemOutput('output')
            with oechem.oemolistream('complex.pdb') as ifs:
                if not oechem.OEReadMolecule(ifs, complex_mol):
                    raise RuntimeError("Error reading complex pdb")
                with oechem.oemolostream('complex.oeb.gz') as ofs:
                    complex_mol.SetData(oechem.OEGetTag('system'), output.encode(system))
                    oechem.OEWriteMolecule(ofs, complex_mol)
            self.success.emit(complex_mol)


        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)

    def end(self):
        #Clean up
        os.remove('protein.pdb')
        os.remove('mol_parameters.ffxml')
        os.remove('pl_tmp.pdb')
        os.remove('nomol_tmp.pdb')
        os.remove('complex.pdb')


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
        default=1000,
        help_text="Number of MD steps")

    complex_mol = parameter.DataSetInputParameter(
        'complex_mol',
        default='complex.oeb.gz',
        #required=True,
        help_text='Single protein to Dock Against')

    def begin(self):
        # Initialize openmm integrator
        self.integrator = openmm.LangevinIntegrator(self.args.temperature*unit.kelvin, 1/unit.picoseconds, 0.002*unit.picoseconds)

    def process(self, mol, port):
        try:
            self.log.info('Regenerating positions and topology from OEMol')
            def extractPositionsFromOEMOL(molecule):
                # Taken from openmoltools
                positions = unit.Quantity(np.zeros([molecule.NumAtoms(), 3], np.float32), unit.angstroms)
                coords = molecule.GetCoords()
                for index in range(molecule.NumAtoms()):
                    positions[index,:] = unit.Quantity(coords[index], unit.angstroms)
                return positions
            positions = extractPositionsFromOEMOL(mol)
            topology = forcefield_generators.generateTopologyFromOEMol(mol)

            # Reconstruct the OpenMM system
            if 'system' in mol.GetData().keys():
                self.log.info('Reconstructing System from mol')
                serialized_system = mol.GetData(oechem.OEGetTag('system'))
                system = openmm.XmlSerializer.deserialize( serialized_system )
                # Initialize Simulation
                simulation = app.Simulation(topology, system, self.integrator)
                #simulation = app.Simulation(topology, system, self.integrator, openmm.Platform.getPlatformByName('CPU'))
                platform = simulation.context.getPlatform().getName()
                self.log.info('Running OpenMMSimulation on Platform {}'.format(platform))
            else:
                raise RuntimeError('Could not find system from mol')

            # Check if mol has State data attached
            if 'state' in mol.GetData().keys():
                outfname = 'restart'
                mol.GetData(oechem.OEGetTag('state'))
                serialized_state = mol.GetData(oechem.OEGetTag('state'))
                state = openmm.XmlSerializer.deserialize( serialized_state )
                simulation.context.setState(state)
            else:
                self.log.info('Minimizing system...')
                outfname = 'simulation'
                # Set initial positions and velocities then minimize
                simulation.context.setPositions(positions)
                simulation.context.setVelocitiesToTemperature(self.args.temperature*unit.kelvin)
                # Temporarily, place some restrictions on minization to run faster
                #simulation.minimizeEnergy()
                simulation.minimizeEnergy(tolerance=unit.Quantity(10.0,unit.kilojoules/unit.moles),maxIterations=20)
                st = simulation.context.getState(getPositions=True,getEnergy=True)
                self.log.info('\tMinimized energy is {}'.format(st.getPotentialEnergy()))

            # Do MD simulation and report energies
            statereporter = app.StateDataReporter(outfname+'.log', 100, step=True, potentialEnergy=True, temperature=True)
            simulation.reporters.append(statereporter)
            outlog = open(outfname+'.log', 'r')
            simulation.step(self.args.steps)
            self.log.info(outlog.read())

            # Save serialized State object
            state = simulation.context.getState( getPositions=True,
                                              getVelocities=True,
                                              getParameters=True )

            # Attach openmm objects to mol, emit to output
            output = OpenMMSystemOutput('output')
            with oechem.oemolostream(outfname+'.oeb.gz') as ofs:
                mol.SetData(oechem.OEGetTag('system'), output.encode(system))
                mol.SetData(oechem.OEGetTag('state'), output.encode(state))
                mol.SetData(oechem.OEGetTag('log'), outlog.read())
                oechem.OEWriteMolecule(ofs, mol)
            self.success.emit(mol)

        except Exception as e:
                # Attach error message to the molecule that failed
                self.log.error(traceback.format_exc())
                mol.SetData('error', str(e))
                # Return failed mol
                self.failure.emit(mol)
