import io, os, traceback
import uuid
import numpy as np
import mdtraj, parmed
from openeye import oechem
from simtk import openmm, unit
from simtk.openmm import app
import OpenMMCubes.simtools as simtools
import OpenMMCubes.utils as utils
from ComplexPrepCubes import utils as complex_utils
from floe.api import ParallelOEMolComputeCube, parameter


class OpenMMComplexSetup(ParallelOEMolComputeCube):
    title = "OpenMM Complex Setup"
    version = "0.0.2"
    classification = [["Protein Preparation", "OpenMM", "Forcefield Assignment"],
    ["Protein Preparation", "PDBFixer", "Solvate"],
    ["Protein Preparation", "PDBFixer", "Add Missing Atoms"],
    ["Protien Preparation", "PDBFixer", "Assign Protonation States"],
    ["Protein-Ligand Preparation", "ParmEd", "Generate Complex"]]
    tags = ['PDBFixer', 'OpenMM', 'ParmEd', 'Parallel Cube']
    description = """
    Using PDBFixer, add missing atoms, assign protonation state with a given pH,
    solvate the system with TIP3P, and assign forcefield parameters (default: amber99sbildn).
    Generate a parameterized parmed Structure of the solvated protein:ligand complex.

    Input:
    -------
    protein - Requires a PDB file of the protein.
    oechem.OEMol - Streamed-in charged and docked molecule with explicit hydrogens.

    Output:
    -------
    oechem.OEMol - Emits molecule with attachments:
        - SDData Tags: { Structure: str <parmed.Structure> }
        - Generic Tags: { Structure : parmed.Structure (base64-encoded) }
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 3600}, # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

    protein = parameter.DataSetInputParameter(
        'protein',
        required=True,
        help_text='Protein PDB file')

    pH = parameter.DecimalParameter(
        'pH',
        default=7.4,
        help_text="Solvent pH used to select appropriate protein protonation state.")

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10,
        help_text="Padding around protein for solvent box (angstroms)")

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=100,
        help_text="Salt concentration (millimolar)")

    protein_forcefield = parameter.DataSetInputParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        help_text='Forcefield parameters for protein')

    solvent_forcefield = parameter.DataSetInputParameter(
        'solvent_forcefield',
        default='tip3p.xml',
        help_text='Forcefield parameters for solvent')

    def begin(self):
        protein = oechem.OEMol()
        self.args.protein = utils.download_dataset_to_file(self.args.protein)
        # Read the PDB file into an OEMol
        with oechem.oemolistream(self.args.protein) as ifs:
            if not oechem.OEReadMolecule(ifs, protein):
                raise RuntimeError("Error reading protein")
        # Read the PDB file into an OpenMM PDBFile object
        self.proteinpdb = app.PDBFile(self.args.protein)

        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, mol, port):
        try:
            # Check for generic data.
            if utils.PackageOEMol.checkTags(mol, ['Structure']):
                gd = utils.PackageOEMol.unpack(mol)
                molecule_structure = gd['Structure']
                self.opt['outfname'] = '{}-complex'.format(gd['IDTag'])

            # Generate parameterized protein Structure
            protein_structure = simtools.genProteinStructure(
                self.proteinpdb, **self.opt)

            # Merge structures to prevent adding solvent in pocket
            # Ligand must be in docked position
            pl_structure = simtools.mergeStructure(
                protein_structure, molecule_structure)
            self.log.info('{}: {}'.format(self.opt['outfname'], pl_structure))

            # Returns solvated system w/o ligand.
            solv_structure = simtools.solvateComplexStructure(
                pl_structure, **self.opt)

            # Remerge with ligand structure
            full_structure = simtools.mergeStructure(
                solv_structure, molecule_structure)
            self.log.info('Solvated {}: {}'.format(
                self.opt['outfname'], full_structure))
            self.log.info('\tBox = {}'.format(full_structure.box))

            
            # Attach the Structure to OEMol
            oechem.OESetSDData(mol, 'Structure', str(full_structure))
            packedmol = utils.PackageOEMol.pack(mol, full_structure)
            packedmol.SetData(oechem.OEGetTag(
                'outfname'), self.opt['outfname'])

            # Save the reference positions in OEMol
            ref_positions = full_structure.positions
            packedpos = utils.PackageOEMol.encodePyObj(ref_positions)
            packedmol.SetData(oechem.OEGetTag('OEMDDataRefPositions'), packedpos)

            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mol)


class OpenMMminimizeCube(ParallelOEMolComputeCube):
    title = 'Minimize the molecule system by using OpenMM'

    version = "0.0.0"

    classification = [["Simulation", "OpenMM", "Minimization"]]
    tags = ['OpenMM', 'MDTraj', 'Parallel Cube']

    description = """
    Minimize the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and minimize it.

    Input parameters:
    steps (integer): the number of steps of minimization to apply. If 0
    the minimization will procede until convergence is reached

    Input parameters:
    """
    
    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 28800}, # Default 8 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

    steps = parameter.IntegerParameter(
        'steps',
        default=0,
        help_text="""Number of minimization steps. 
                  If 0 the mimimization will continue 
                  until convergence""")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text="Mask selection to apply restraints")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=5.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol ang^2)")

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300,
        help_text="Temperature (Kelvin)")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald.")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The nonbonded cutoff in angstroms.
        This is ignored if nonbondedMethod is NoCutoff.""")

    constraints = parameter.StringParameter(
        'constraints',
        default='HBonds',
        choices=['None', 'HBonds', 'HAngles', 'AllBonds'],
        help_text="""None, HBonds, HAngles, or AllBonds
        Which type of constraints to add to the system (e.g., SHAKE).
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained""")

    outfname = parameter.StringParameter(
        'outfname',
        default='min',
        help_text='Filename suffix for output simulation files. Formatted: <title>-<outfname>')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        description='Increase log file verbosity')

    platform =  parameter.StringParameter(
        'platform',
        default='Auto', 
        choices=['Auto', 'Reference', 'CPU', 'CUDA', 'OpenCL'],
        help_text='Select which platform to use to run the simulation')

    
    def begin(self):
        self.opt = vars( self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'min'

        return

    def process(self, mol, port):
        try:
            if utils.PackageOEMol.checkTags(mol, ['Structure']):
                gd = utils.PackageOEMol.unpack(mol)
                self.opt['outfname'] = '{}-{}'.format(gd['IDTag'], self.opt['outfname'])

            mdData = utils.MDData(mol)

            self.opt['molecule'] = mol

            self.log.info('MINIMIZING System: %s' % gd['IDTag'])
            simtools.simulation(mdData, **self.opt)
        
            packedmol = mdData.packMDData(mol)

            # Update the OEMol complex positions to match the new
            # Parmed structure
            new_temp_mol = complex_utils.openmmTop_to_oemol(mdData.topology, mdData.positions)
            new_pos = new_temp_mol.GetCoords()
            packedmol.SetCoords(new_pos)

            # packedmol.SetData(oechem.OEGetTag(
            #     'outfname'), self.opt['outfname'])

            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)
   
        return

            
class OpenMMnvtCube(ParallelOEMolComputeCube):
    title = 'OpenMM NVT simulation'
    version = "0.0.0"
    classification = [["Simulation", "OpenMM", "NVT"]]
    tags = ['OpenMM', 'MDTraj', 'Parallel Cube']

    description = """NVT simulation of the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and will perform a MD simulation
    at constant temperature and volume

    Input parameters:
    ----------------
      picosec (decimal): Number of picoseconds to warm up the complex.
      temperature (decimal): target temperature.
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 28800}, # Default 8 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }
    
    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    time = parameter.DecimalParameter(
        'time',
        default=10.0,
        help_text="NVT simulation time in picoseconds")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text="Mask selection to apply restraints")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol ang^2)")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The nonbonded cutoff in angstroms.
        This is ignored if nonbondedMethod is NoCutoff.""")

    constraints = parameter.StringParameter(
        'constraints',
        default='HBonds',
        choices=['None', 'HBonds', 'HAngles', 'AllBonds'],
        help_text="""None, HBonds, HAngles, or AllBonds
        Which type of constraints to add to the system (e.g., SHAKE).
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained""")

    trajectory_filetype = parameter.StringParameter(
        'trajectory_filetype',
        default='NetCDF',
        choices=['NetCDF', 'DCD', 'HDF5'],
        help_text="NetCDF, DCD, HDF5. Filetype to write trajectory files")

    trajectory_selection = parameter.StringParameter(
        'trajectory_selection',
        default=None,
        choices=[None, 'protein or resname LIG', 'protein', 'resname LIG'],
        help_text='atoms subset to write in trajectory')

    trajectory_interval = parameter.IntegerParameter(
        'trajectory_interval',
        default=1000,
        help_text="Step interval for trajetory snapshots.")

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=10000,
        help_text="Step interval for reporting data.")

    outfname = parameter.StringParameter(
        'outfname',
        default='nvt',
        help_text='Filename suffix for output simulation files. Formatted: <title>-<outfname>')

    tarxz = parameter.BooleanParameter(
        'tarxz',
        default=True,
        description='Create a tar.xz file of the attached data')

    center = parameter.BooleanParameter(
        'center',
        default=True,
        description='Center the system to the OpenMM unit cell')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        description='Increase log file verbosity.')

    platform = parameter.StringParameter(
        'platform',
        default='Auto',
        choices=['Auto', 'Reference', 'CPU', 'CUDA', 'OpenCL'],
        help_text='Select which platform to use to run the simulation')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['convert'] = False
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'nvt'
        
        conv_rule = [self.opt['trajectory_selection'] != None,
                     self.opt['trajectory_filetype'] != 'NetCDF']
        
        if any(conv_rule):
            self.opt['convert'] = True

        return

    def process(self, mol, port):
        try:
            if utils.PackageOEMol.checkTags(mol, ['Structure']):
                gd = utils.PackageOEMol.unpack(mol)
                self.opt['outfname'] = '{}-{}'.format(gd['IDTag'], self.opt['outfname'])

            mdData = utils.MDData(mol)

            self.opt['molecule'] = mol

            self.log.info('START NVT SIMULATION: %s' % gd['IDTag'])
            simtools.simulation(mdData, **self.opt)
                
            packedmol = mdData.packMDData(mol)

            # Update the OEMol complex positions to match the new
            # Parmed structure
            new_temp_mol = complex_utils.openmmTop_to_oemol(mdData.topology, mdData.positions)
            new_pos = new_temp_mol.GetCoords()
            packedmol.SetCoords(new_pos)

            # Create a tar.xz archive of the generic data and trajectories
            if self.opt['tarxz']:
                utils.PackageOEMol.dump(
                    packedmol, outfname=self.opt['outfname'], tarxz=self.opt['tarxz'])
            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)
         
        return

    
class OpenMMnptCube(ParallelOEMolComputeCube):
    title = 'OpenMM NPT simulation'
    version = "0.0.0"
    classification = [["Simulation", "OpenMM", "NPT"]]
    tags = ['OpenMM', 'MDTraj', 'Parallel Cube']

    description = """NPT simulation of the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and will perform a MD simaultion at
    constant tempertaure and pressure.

    Input parameters:
    ----------------
      picosec (decimal): Number of picoseconds to perform the complex simulation.
      temperature (decimal): target temperature.
      pressure (decimal): target pressure.
    """
    
    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 28800}, # Default 8 hour limit (units are seconds)
        "item_count": {"default": 1} # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default= 300,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default= 1.0,
        help_text="Pressure (atm)")

    time = parameter.DecimalParameter(
        'time',
        default=10.0,
        help_text="NPT simulation time in picoseconds")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text="Mask selection to apply restraints")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol ang^2)")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald.")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The nonbonded cutoff in angstroms.
        This is ignored if nonbondedMethod is NoCutoff.""")

    constraints = parameter.StringParameter(
        'constraints',
        default='HBonds',
        choices=['None', 'HBonds', 'HAngles', 'AllBonds'],
        help_text="""None, HBonds, HAngles, or AllBonds
        Which type of constraints to add to the system (e.g., SHAKE).
        None means no bonds are constrained.
        HBonds means bonds with hydrogen are constrained""")

    trajectory_filetype = parameter.StringParameter(
        'trajectory_filetype',
        default='NetCDF',
        choices=['NetCDF', 'DCD', 'HDF5'],
        help_text="NetCDF, DCD, HDF5. Filetype to write trajectory files")

    trajectory_selection = parameter.StringParameter(
        'trajectory_selection',
        default=None,
        choices=[None, 'protein or resname LIG', 'protein', 'resname LIG'],
        help_text='atoms subset to write in trajectory')

    trajectory_interval = parameter.IntegerParameter(
        'trajectory_interval',
        default=1000,
        help_text="Step interval for trajetory snapshots.")

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=10000,
        help_text="Step interval for reporting data.")

    outfname = parameter.StringParameter(
        'outfname',
        default='npt',
        help_text='Filename suffix for output simulation files. Formatted: <title>-<outfname>')

    tarxz = parameter.BooleanParameter(
        'tarxz',
        default=True,
        description='Create a tar.xz file of the attached data')

    center = parameter.BooleanParameter(
        'center',
        default=True,
        description='Center the system to the OpenMM unit cell')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        description='Increase log file verbosity.')

    platform = parameter.StringParameter(
        'platform',
        default='Auto', 
        choices=[ 'Auto', 'Reference', 'CPU', 'CUDA', 'OpenCL'],
        help_text= 'Select which platform to use to run the simulation')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['convert'] = False
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'npt'
        conv_rule = [self.opt['trajectory_selection'] != None,
                     self.opt['trajectory_filetype'] != 'NetCDF']
        
        if any(conv_rule):
            self.opt['convert'] = True

        return

    def process(self, mol, port):
        try:
            if utils.PackageOEMol.checkTags(mol, ['Structure']):
                gd = utils.PackageOEMol.unpack(mol)
                self.opt['outfname'] = '{}-{}'.format(gd['IDTag'], self.opt['outfname'])

            mdData = utils.MDData(mol)

            self.opt['molecule'] = mol

            self.log.info('START NPT SIMULATION %s' % gd['IDTag'])
            simtools.simulation(mdData, **self.opt)
                
            packedmol = mdData.packMDData(mol)

            # Update the OEMol complex positions to match the new
            # Parmed structure
            new_temp_mol = complex_utils.openmmTop_to_oemol(mdData.topology, mdData.positions)
            new_pos = new_temp_mol.GetCoords()
            packedmol.SetCoords(new_pos)

            # Create a tar.xz archive of the generic data and trajectories
            if self.opt['tarxz']:
                utils.PackageOEMol.dump(
                    packedmol, outfname=self.opt['outfname'], tarxz=self.opt['tarxz'])
            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)
        return
