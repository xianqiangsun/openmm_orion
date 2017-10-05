import traceback
import OpenMMCubes.simtools as simtools
import OpenMMCubes.utils as utils
from floe.api import ParallelOEMolComputeCube, parameter
from openeye import oechem


class OpenMMminimizeCube(ParallelOEMolComputeCube):
    title = 'Minimization Cube'

    version = "0.0.0"
    classification = [["Simulation", "OpenMM", "Minimization"]]
    tags = ['OpenMM', 'Parallel Cube']

    description = """
    Minimize the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and minimize it.

    Input parameters:
    steps (integer): the number of steps of minimization to apply. If 0
    the minimization will proceed until convergence is reached
    """
    
    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    steps = parameter.IntegerParameter(
        'steps',
        default=0,
        help_text="""Number of minimization steps. 
                  If 0 the minimization will continue 
                  until convergence""")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text="""Mask selection to apply restraints. Possible keywords are:
                  ligand, protein, water, ions, ca_protein, cofactors. 
                  The selection can be refined by using logical tokens: 
                  not, noh, and, or, diff, around""")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=5.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)")

    freeze = parameter.StringParameter(
        'freeze',
        default='',
        help_text="""Mask selection to freeze atoms along the MD 
                  simulation. Possible keywords are: ligand, protein, water, 
                  ions, ca_protein, cofactors. The selection can be refined by
                  using logical tokens: not, noh, and, or, diff, around""")

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300,
        help_text="Temperature (Kelvin)")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if the non-bonded method is NoCutoff.""")

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
        help_text='Filename suffix for output simulation files')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        description='Increase log file verbosity')

    platform = parameter.StringParameter(
        'platform',
        default='Auto', 
        choices=['Auto', 'Reference', 'CPU', 'CUDA', 'OpenCL'],
        help_text='Select which platform to use to run the simulation')

    cuda_opencl_precision = parameter.StringParameter(
        'cuda_opencl_precision',
        default='single',
        choices=['single', 'mixed', 'double'],
        help_text='Select the CUDA or OpenCL precision')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'min'

        return

    def process(self, mol, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(mol) if dp.GetTag() in
                        ["temperature"]}
            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(mol.GetTitle(), new_args))
                opt.update(new_args)

            if utils.PackageOEMol.checkTags(mol, ['Structure']):
                gd = utils.PackageOEMol.unpack(mol)
                opt['outfname'] = '{}-{}'.format(gd['IDTag'], self.opt['outfname'])

            mdData = utils.MDData(mol)

            opt['molecule'] = mol

            self.log.info('MINIMIZING System: %s' % gd['IDTag'])
            simtools.simulation(mdData, **opt)

            packedmol = mdData.packMDData(mol)

            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)
   
        return

            
class OpenMMnvtCube(ParallelOEMolComputeCube):
    title = 'NVT Cube'
    version = "0.0.0"
    classification = [["Simulation", "OpenMM", "NVT"]]
    tags = ['OpenMM', 'Parallel Cube']

    description = """NVT simulation of the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and will perform a MD simulation
    at constant temperature and volume

    Input parameters:
    ----------------
      picosec (decimal): Number of picoseconds to warm up the complex.
      temperature (decimal): target temperature
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
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
        help_text=""""Mask selection to apply restraints. Possible keywords are: 
                  ligand, protein, water, ions, ca_protein, cofactors. 
                  The selection can be refined by using logical tokens: 
                  not, noh, and, or, diff, around""")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol A^2)")

    freeze = parameter.StringParameter(
        'freeze',
        default='',
        help_text="""Mask selection to freeze atoms along the MD 
                  simulation. Possible keywords are: ligand, protein, water, 
                  ions, ca_protein, cofactors. The selection can be refined by
                  using logical tokens: not, noh, and, or, diff, around""")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if non-bonded method is NoCutoff.""")

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
        default='DCD',
        choices=['DCD', 'NetCDF', 'HDF5'],
        help_text="NetCDF, DCD, HDF5. File type to write trajectory files")

    trajectory_interval = parameter.IntegerParameter(
        'trajectory_interval',
        default=0,
        help_text="Step interval for trajectory snapshots. If 0 the trajectory"
                  "file will not be generated")

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=0,
        help_text="Step interval for reporting data. If 0 the reporter file"
                  "will not be generated")

    outfname = parameter.StringParameter(
        'outfname',
        default='nvt',
        help_text='Filename suffix for output simulation files')

    tar = parameter.BooleanParameter(
        'tar',
        default=False,
        description='Create a tar.xz file of the attached data')

    center = parameter.BooleanParameter(
        'center',
        default=False,
        description='Center the system to the OpenMM unit cell')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=True,
        description='Increase log file verbosity')

    platform = parameter.StringParameter(
        'platform',
        default='Auto',
        choices=['Auto', 'Reference', 'CPU', 'CUDA', 'OpenCL'],
        help_text='Select which platform to use to run the simulation')

    cuda_opencl_precision = parameter.StringParameter(
        'cuda_opencl_precision',
        default='single',
        choices=['single', 'mixed', 'double'],
        help_text='Select the CUDA or OpenCL precision')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'nvt'

        return

    def process(self, mol, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(mol) if dp.GetTag() in
                        ["temperature"]}
            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(mol.GetTitle(), new_args))
                opt.update(new_args)

            if utils.PackageOEMol.checkTags(mol, ['Structure']):
                gd = utils.PackageOEMol.unpack(mol)
                opt['outfname'] = '{}-{}'.format(gd['IDTag'], self.opt['outfname'])

            mdData = utils.MDData(mol)

            opt['molecule'] = mol

            self.log.info('START NVT SIMULATION: %s' % gd['IDTag'])

            simtools.simulation(mdData, **opt)
                
            packedmol = mdData.packMDData(mol)

            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)
         
        return

    
class OpenMMnptCube(ParallelOEMolComputeCube):
    title = 'NPT Cube'
    version = "0.0.0"
    classification = [["Simulation", "OpenMM", "NPT"]]
    tags = ['OpenMM', 'Parallel Cube']

    description = """NPT simulation of the protein:ligand complex.

    This cube will take in the streamed complex.oeb.gz file containing
    the solvated protein:ligand complex and will perform a MD simulation at
    constant temperature and pressure.

    Input parameters:
    ----------------
      picosec (decimal): Number of picoseconds to perform the complex simulation.
      temperature (decimal): target temperature
      pressure (decimal): target pressure
    """
    
    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1}, # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    time = parameter.DecimalParameter(
        'time',
        default=10.0,
        help_text="NPT simulation time in picoseconds")

    restraints = parameter.StringParameter(
        'restraints',
        default='',
        help_text=""""Mask selection to apply restraints. Possible keywords are: 
                  ligand, protein, water, ions, ca_protein, cofactors. 
                  The selection can be refined by using logical tokens: 
                  not, noh, and, or, diff, around""")

    restraintWt = parameter.DecimalParameter(
        'restraintWt',
        default=2.0,
        help_text="Restraint weight for xyz atom restraints in kcal/(mol ang^2)")

    freeze = parameter.StringParameter(
        'freeze',
        default='',
        help_text="""Mask selection to freeze atoms along the MD simulation. 
                  Possible keywords are: ligand, protein, water, ions, ca_protein, 
                  cofactors. The selection can be refined by using logical tokens: 
                  not, noh, and, or, diff, around""")

    nonbondedMethod = parameter.StringParameter(
        'nonbondedMethod',
        default='PME',
        choices=['NoCutoff', 'CutoffNonPeriodic',
                 'CutoffPeriodic', 'PME', 'Ewald'],
        help_text="NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald.")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10,
        help_text="""The non-bonded cutoff in angstroms.
        This is ignored if non-bonded method is NoCutoff""")

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
        default='DCD',
        choices=['DCD', 'NetCDF', 'HDF5'],
        help_text="NetCDF, DCD, HDF5. File type to write trajectory files")

    trajectory_interval = parameter.IntegerParameter(
        'trajectory_interval',
        default=0,
        help_text="Step interval for trajectory snapshots. If 0 the trajectory"
                  "file will not be generated")

    reporter_interval = parameter.IntegerParameter(
        'reporter_interval',
        default=0,
        help_text="Step interval for reporting data. If 0 the reporter file"
                  "will not be generated")

    outfname = parameter.StringParameter(
        'outfname',
        default='npt',
        help_text='Filename suffix for output simulation files. Formatted: <title>-<outfname>')

    tar = parameter.BooleanParameter(
        'tar',
        default=False,
        description='Create a tar.xz file of the attached data')

    center = parameter.BooleanParameter(
        'center',
        default=False,
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

    cuda_opencl_precision = parameter.StringParameter(
        'cuda_opencl_precision',
        default='single',
        choices=['single', 'mixed', 'double'],
        help_text='Select the CUDA or OpenCL precision')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.opt['SimType'] = 'npt'

        return

    def process(self, mol, port):
        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(mol) if dp.GetTag() in
                        ["temperature", "pressure"]}

            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(mol.GetTitle(), new_args))

                opt.update(new_args)

            if utils.PackageOEMol.checkTags(mol, ['Structure']):
                gd = utils.PackageOEMol.unpack(mol)
                opt['outfname'] = '{}-{}'.format(gd['IDTag'], self.opt['outfname'])

            mdData = utils.MDData(mol)

            opt['molecule'] = mol

            self.log.info('START NPT SIMULATION %s' % gd['IDTag'])
            simtools.simulation(mdData, **opt)
                
            packedmol = mdData.packMDData(mol)

            self.success.emit(packedmol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mol.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(mol)

        return