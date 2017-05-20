import io, os, sys, base64, parmed, mdtraj, pdbfixer
import numpy as np
from sys import stdout
from openeye import oechem
from simtk import unit, openmm
from simtk.openmm import app
import OpenMMCubes.utils as utils
import pyparsing as pyp


try:
    import cPickle as pickle
except ImportError:
    import pickle

    
def genProteinStructure(proteinpdb, **opt):
    """
    Starting from OpenMM PDBFile, generates the OpenMM System for the protein and
    then create the parametrized ParmEd Structure of the protein.

    Parameters
    ----------
    proteinpdb : openmm.app.PDBFile object,
        Loaded PDBFile object of the protein.
    protein_forcefield : (opt), xml file, default='amber99sbildn.xml'
        Forcefield parameters for protein
    solvent_forcefield : opt), xml file, default='tip3p.xml'
        Forcefield parameters for solvent

    Returns
    -------
    solv_structure : parmed.structure.Structure
        The parameterized Structure of the protein with solvent molecules. (No ligand).
    """
    #Generate protein Structure object
    forcefield = app.ForceField(opt['protein_forcefield'], opt['solvent_forcefield'])
    protein_system = forcefield.createSystem( proteinpdb.topology )
    protein_structure = parmed.openmm.load_topology(proteinpdb.topology,
                                                    protein_system,
                                                    xyz=proteinpdb.positions)
    return protein_structure

def combinePositions(proteinPositions, molPositions):
    """
    Loops through the positions from the ParmEd structures of the protein and ligand,
    divides by unit.angstroms which will ensure both positions arrays are in the same units.

    Parameters
    ----------
    proteinPositions : list of 3-element Quantity tuples.
        Positions list taken directly from the protein Structure.
    molPositions : list of 3-element Quantity tuples.
        Positions list taken directly from the molecule Structure.

    Returns
    -------
    positions : list of 3-element Quantity tuples.
        ex. unit.Quantity(positions, positions_unit)
        Combined positions of the protein and molecule Structures.
    """
    positions_unit = unit.angstroms
    positions0_dimensionless = np.array(proteinPositions / positions_unit)
    positions1_dimensionless = np.array(molPositions / positions_unit)
    coordinates = np.vstack(
        (positions0_dimensionless, positions1_dimensionless))
    natoms = len(coordinates)
    positions = np.zeros([natoms, 3], np.float32)
    for index in range(natoms):
            (x, y, z) = coordinates[index]
            positions[index, 0] = x
            positions[index, 1] = y
            positions[index, 2] = z
    positions = unit.Quantity(positions, positions_unit)
    return positions

def mergeStructure(proteinStructure, molStructure):
    """
    Combines the parametrized ParmEd structures of the protein and ligand to
    create the Structure for the protein:ligand complex, while retaining the SMIRFF
    parameters on the ligand. Preserves positions and box vectors.
    (Not as easily achieved using native OpenMM tools).

    Parameters
    ----------
    proteinStructure : parmed.structure.Structure
        The parametrized structure of the protein.
    moleculeStructure : parmed.structure.Structure
        The parametrized structure of the ligand.

    Returns
    -------
    structure : parmed.structure.Structure
        The parametrized structure of the protein:ligand complex.
    """
    structure = proteinStructure + molStructure
    positions = combinePositions(proteinStructure.positions, molStructure.positions)
    # Concatenate positions arrays (ensures same units)
    structure.positions = positions
    # Restore original box vectors
    structure.box = proteinStructure.box
    return structure

def solvateComplexStructure(structure, **opt):
    """
    Uses PDBFixer to add missing atoms, assign protonation states, and solvate
    the parametrized ParmEd Structure of the protein:ligand complex. Returns the
    parametrized Structure without the ligand present.
    We initially combine the ligand with the protein to prevent PDBFixer
    from adding solvent into the binding pocket.

    Parameters
    ----------
    structure : parmed.structure.Structure
        The parametrized structure of the protein:ligand complex
    pH : (opt), float, default=7.4
        Solvent pH used to select appropriate protein protonation state.
    solvent_padding : (opt), float, default=10
        Padding around protein for solvent box (unit.angstroms)
    salt_concentration : (opt), float, default=10
        Salt concentration (millimolar)
    protein_forcefield : (opt), xml file, default='amber99sbildn9sbildn.xml'
        Forcefield parameters for protein
    solvent_forcefield : opt), xml file, default='tip3p.xml'
        Forcefield parameters for solvent

    Returns
    -------
    solv_structure : parmed.structure.Structure
        The parameterized Structure of the protein with solvent molecules. (No ligand).
    """
    if opt['Logger'] is None:
        printfile = sys.stdout
    else:
        printfile = opt['Logger'].file

    tmpfile = opt['outfname']+'-pl.tmp'
    structure.save(tmpfile,format='pdb')

    seqres = False
    with open(tmpfile, 'r') as infile:
        for line in infile:
            if 'SEQRES' in line:
                seqres = True
                break
    if not seqres:
        print('Did not find SEQRES in PDB. PDBFixer will not find missing Residues.', file=printfile)

    # Solvate with PDBFixer
    print('PDBFixer settings for {}'.format(opt['outfname']), file=printfile)
    print('\tpH = {}'.format(opt['pH']), file=printfile)
    print('\tpadding = {}'.format(unit.Quantity(opt['solvent_padding'], unit.angstroms)), file=printfile)
    print('\tsalt conc. = {}'.format(unit.Quantity(opt['salt_concentration'], unit.millimolar)), file=printfile)
    fixer = pdbfixer.PDBFixer(tmpfile)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()

    if fixer.missingAtoms:
        print('Found missing Atoms:', fixer.missingAtoms, file=printfile)
    if fixer.missingTerminals:
        print('Found missing Terminals:', fixer.missingTerminals, file=printfile)
    if fixer.nonstandardResidues:
        print('Found nonstandard Residues:', fixer.nonstandardResidues, file=printfile)

    fixer.replaceNonstandardResidues()
    #fixer.removeHeterogens(False)
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(opt['pH'])
    fixer.addSolvent(padding=unit.Quantity(opt['solvent_padding'], unit.angstroms),
                ionicStrength=unit.Quantity(opt['salt_concentration'], unit.millimolar))

    # Load PDBFixer object back to Structure
    tmp = parmed.openmm.load_topology(fixer.topology, xyz=fixer.positions)

    # Remove ligand from protein Structure by AmberMask selection
    tmp.strip(":LIG")
    tmp.save(opt['outfname']+'-nomol.tmp',format='pdb')
    # Reload PDBFile
    nomol = app.PDBFile(opt['outfname']+'-nomol.tmp')
    forcefield = app.ForceField(opt['protein_forcefield'], opt['solvent_forcefield'])
    nomol_system = forcefield.createSystem(nomol.topology, rigidWater=False)
    # Regenerate parameterized solvated protein structure
    solv_structure = parmed.openmm.load_topology(nomol.topology,
                                                nomol_system,
                                                xyz=nomol.positions)
    # Restore box vectors
    solv_structure.box = tmp.box

    tmpfiles = [ opt['outfname']+'-pl.tmp', opt['outfname']+'-nomol.tmp' ]
    utils.cleanup(tmpfiles)

    return solv_structure


def simulation(mdData, **opt):
    """
    Minimization, NVT and NPT MD run

    Parameters
    ----------
    mdData : MDData data object
        The object which recovers the Parmed structure data relevant for MD
    """
    
    if opt['Logger'] is None:
        printfile = sys.stdout
    else:
        printfile = opt['Logger'].file

    structure = mdData.structure
    topology = mdData.topology
    positions = mdData.positions
    velocities = mdData.velocities
    box = mdData.box
    stepLen = 0.002

    system = structure.createSystem(nonbondedMethod=eval("app.%s" % opt['nonbondedMethod']),
                                    nonbondedCutoff=opt['nonbondedCutoff']*unit.angstroms,
                                    constraints=eval("app.%s" % opt['constraints']))
    
    integrator = openmm.LangevinIntegrator(opt['temperature']*unit.kelvin, 1/unit.picoseconds, stepLen*unit.picoseconds)
    
    if opt['SimType'] == 'npt':
        # Add Force Barostat to the system
        system.addForce(openmm.MonteCarloBarostat(opt['pressure']*unit.atmospheres, opt['temperature']*unit.kelvin, 25))

    if opt['restraints']:
        opt['Logger'].info("RESTRAINTS mask applied to: {}".format(opt['restraints']))
        #Select atom to restrain
        res_atom_set = restraints(structure, res_mask=opt['restraints'])
        opt['Logger'].info("Number of restainst atoms: {}".format(len(res_atom_set)))
        # define the custom force to restrain atoms to their starting positions
        force_restr = openmm.CustomExternalForce('k_restr*( (x-x0)^2 + (y-y0)^2 + (z-z0)^2 )')
        # Add the restraint weight as a global parameter in kcal/mol/A^2
        force_restr.addGlobalParameter("k_restr", opt['restraintWt']*unit.kilocalories_per_mole/unit.angstroms**2)
        # Define the target xyz coords for the restraint as per-atom (per-particle) parameters
        force_restr.addPerParticleParameter("x0")
        force_restr.addPerParticleParameter("y0")
        force_restr.addPerParticleParameter("z0")

        for idx in range(0,len(positions)):
            if idx in res_atom_set:
                xyz = positions[idx].in_units_of(unit.nanometers)/unit.nanometers
                force_restr.addParticle(idx, xyz)
        
        system.addForce(force_restr)

        
    if opt['platform'] == 'Auto':
        simulation = app.Simulation(topology, system, integrator)
    else:
        try :
            platform = openmm.Platform.getPlatformByName(opt['platform'])
            simulation = app.Simulation(topology, system, integrator, platform)
        except Exception as e:
            raise ValueError('The selected platform is invalid: {}'.format(str(e)))
    
    # Set starting position and velocities
    simulation.context.setPositions(positions)

    # Set Box dimension
    simulation.context.setPeriodicBoxVectors(box[0],box[1],box[2])

    if opt['SimType'] in ['nvt', 'npt']:
        if velocities is not None:
            opt['Logger'].info('RESTARTING simulaiton from the previous State')
            simulation.context.setVelocities(velocities)
        else:
            # Set the velocities drawing from the Boltzmann distribution at the selected temperature
            opt['Logger'].info('GENERATING a new starting State')
            simulation.context.setVelocitiesToTemperature(opt['temperature']*unit.kelvin)

        # Convert simulation time in steps
        opt['steps'] = int(round(opt['time']/stepLen))
        
        # Set Reportes
        for rep in getReporters(**opt):
            simulation.reporters.append(rep)
            
    # OpenMM platform information
    mmver = openmm.version.version
    mmplat = simulation.context.getPlatform()

    if opt['verbose']:
        # Host information
        from platform import uname
        for k,v in uname()._asdict().items():
            print(k, ':', v, file=printfile)
        # Platform properties
        for prop in mmplat.getPropertyNames():
            val = mmplat.getPropertyValue(simulation.context, prop)
            print(prop, ':', val, file=printfile)
            
    print('OpenMM({}) simulation generated for {} platform'.format(mmver, mmplat.getName()), file=printfile)

    if opt['SimType'] in ['nvt', 'npt']:

        opt['Logger'].info('Running {time} ps = {steps} steps of {SimType} MD at {temperature}K'.format( **opt))
        
        # Start Simulation
        simulation.step(opt['steps'])

        if opt['convert']:
            opt['Logger'].info('Converting trajectories to: {trajectory_filetype}'.format(**opt))
            simtools.mdTrajConvert(simulation, outfname=opt['outfname'],
                               trajectory_selection=opt['trajectory_selection'],
                               trajectory_filetype=opt['trajectory_filetype'])

        state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True)
        
    elif opt['SimType'] == 'min':
        # Start Simulation
        
        state = simulation.context.getState(getEnergy=True)

        print('Initial energy = {}'.format(state.getPotentialEnergy()), file=printfile)

        simulation.minimizeEnergy(maxIterations=opt['steps'])

        state = simulation.context.getState(getPositions=True, getEnergy=True)
        print('Minimized energy = {}'.format(state.getPotentialEnergy()), file=printfile)

    
    # OpenMM Quantity object
    structure.positions = state.getPositions(asNumpy=False)
    # OpenMM Quantity object
    structure.box_vectors = state.getPeriodicBoxVectors()

    if opt['SimType'] in ['nvt', 'npt']:
        # numpy array in units of angstrom/picosecond
        structure.velocities = state.getVelocities(asNumpy=False)
 
    return



def getReporters(totalSteps=None, outfname=None, **opt):
    """
    Creates 3 OpenMM Reporters for the simulation.

    Parameters
    ----------
    totalSteps : int
        The total number of simulation steps
    reportInterval : (opt), int, default=1000
        Step frequency to write to reporter file.
    outfname : str
        Specifies the filename prefix for the reporters.

    Returns
    -------
    reporters : list of three openmm.app.simulation.reporters
        (0) state_reporter: writes energies to '.log' file.
        (1) progress_reporter: prints simulation progress to 'sys.stdout'
        (2) traj_reporter: writes trajectory to '.nc' file. AMBER NetCDF(3.0)
    """
    if totalSteps is None:
        totalSteps = opt['steps']
    if outfname is None:
        outfname = opt['outfname']
    progress_reporter = app.StateDataReporter(stdout, separator="\t",
                                        reportInterval=opt['reporter_interval'],
                                        step=True, totalSteps=totalSteps,
                                        time=True, speed=True, progress=True,
                                        elapsedTime=True, remainingTime=True)

    state_reporter = app.StateDataReporter(outfname+'.log', separator="\t",
                                        reportInterval=opt['reporter_interval'],
                                        step=True,
                                        potentialEnergy=True, totalEnergy=True,
                                        volume=True, temperature=True)

    # Uses NetCDF(4.0) but not VMD compatible.
    #traj_reporter = mdtraj.reporters.HDF5Reporter(outfname+'.h5', opt['reporter_interval'])

    # Default to NetCDF since VMD compatible.
    traj_reporter = mdtraj.reporters.NetCDFReporter(outfname+'.nc', opt['trajectory_interval'])

    reporters = [state_reporter, progress_reporter, traj_reporter]
    return reporters


def mdTrajConvert(simulation=None, outfname=None, trajectory_selection=None, trajectory_filetype='NetCDF'):
    """
    Used to convert the (mdTraj) trajectory file between: HDF5, DCD, or NetCDF.
    Can also be used to write out a subset of the trajectory.
    See mdTraj atom selection docs:
        http://mdtraj.org/1.8.0/examples/atom-selection.html

    Parameters
    ----------
    simulation : openmm.app.simulation.Simulation
        The OpenMM Simulation object used to generate the mdTraj topology for
        writing to DCD or NetCDF files.
    trajectory_selection : str, default=None (writes all atoms)
        Examples: 'protein', 'resname LIG', 'protein or resname LIG'
        See mdTraj atom selection docs:
            http://mdtraj.org/1.8.0/examples/atom-selection.html
    trajectory_filetype : str, default='NetCDF'
        Can be HDF5, DCD, or NetCDF.
    outfname : str
        Specifies the filename prefix for the trajectory files.

    """
    if opt['Logger'] is None:
        printfile = sys.stdout
    else:
        printfile = opt['Logger'].file

    atom_indices = None
    traj_dict = { 'HDF5' : outfname +'.h5',
                  'DCD' : outfname +'.dcd',
                  'NetCDF' : outfname +'.nc' }
    for rep in simulation.reporters:
        try:
            trajfile = rep._traj_file
            trajfname = trajfile._handle.filename
            if trajectory_selection is not None:
                atom_indices = trajfile.topology.select('%s' % trajectory_selection)
            trajfile.close()
        except Exception as e:
            pass
    top = mdtraj.Topology.from_openmm(simulation.topology)
    outfile = traj_dict.get(trajectory_filetype)
    traj = mdtraj.load(trajfname, atom_indices=atom_indices, top=top)
    print("\tTrajectory subset: '{}'\n\t{}".format(trajectory_selection, traj), file=printfile)
    print("\tConverted trajectory to %s" % (outfile), file=printfile)
    traj.save(outfile)


def restraints(structure, res_mask=''):
    """
    This functions select the atom indexes to apply harmonic restarins
    
    Parameters
    ----------
    structure : Parmed structure object
        The Parmed structure 
     
    res_mask : python string
        A string used to select the atom to restraints. A BNF grammar is defined
        by using the tokens: "ligand", "protein", "water", "ions", "excipients".
        Operator tokens are "not" and "and".
      
    Returns
    -------
    atom_set : python set
        the select atom indexes to apply harmonic restraints
   
    Notes
    -----
    Example of selection string:
        res_mask = "ligand and protein"
        res_mask = "not water and not ions"
        res_mask = "ligand and protein and excipients"
    """

    def build_set(ls, dsets):
        # Terminal Literal return the related set
        if isinstance(ls, str):
            return dsets[ls]
        # Not
        if len(ls) == 2:
            return system - build_set(ls[1], dsets)
        # And
        if len(ls) == 3:
            return build_set(ls[0], dsets) | build_set(ls[2], dsets)
        else:
            raise ValueError("The passed list {} cannot have more than 3 tokens".format(ls))

    # Parse Action-Maker
    def makeLRlike(numterms):
        if numterms is None:
            # None operator can only by binary op
            initlen = 2
            incr = 1
        else:
            initlen = {0: 1, 1: 2, 2: 3, 3: 5}[numterms]
            incr = {0: 1, 1: 1, 2: 2, 3: 4}[numterms]

        # Define parse action for this number of terms,
        # to convert flat list of tokens into nested list
        def pa(s, l, t):
            t = t[0]
            if len(t) > initlen:
                ret = pyp.ParseResults(t[:initlen])
                i = initlen
                while i < len(t):
                    ret = pyp.ParseResults([ret] + t[i:i + incr])
                    i += incr
                return pyp.ParseResults([ret])

        return pa

    operand = pyp.Literal("protein") | pyp.Literal("ligand") | \
              pyp.Literal("water") | pyp.Literal("ions") | \
              pyp.Literal("excipients")

    # BNF Grammar definition with parseAction makeLRlike
    expr = pyp.operatorPrecedence(operand,
                              [
                                  (None, 2, pyp.opAssoc.LEFT, makeLRlike(None)),
                                  (pyp.Literal("not"), 1, pyp.opAssoc.RIGHT, makeLRlike(1)),
                                  (pyp.Literal("and"), 2, pyp.opAssoc.LEFT, makeLRlike(2)),
                              ])

    proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO',
                       'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS',
                       'PHE', 'SER', 'TRP', 'VAL']

    ligand = set()
    protein = set()
    water = set()
    ions = set()
    excipients = set()
    system = set()
    
    for res in structure.residues:
        if res.name == 'LIG':
            for at in res.atoms:
                ligand.add(at.idx)
        elif res.name == 'HOH':
            for at in res.atoms:
                water.add(at.idx)
        elif len(res.atoms) == 1:
            for at in res.atoms:
                ions.add(at.idx)
        elif res.name in proteinResidues:
            for at in res.atoms:
                protein.add(at.idx)
        else:
            for at in res.atoms:
                excipients.add(at.idx)

    dic_set = {'ligand':ligand, 'protein':protein, 'water':water, 'ions':ions, 'excipients':excipients}

    for k in dic_set:
        system = system | dic_set[k]

    # print("ligand = {}".format(len(ligand)))
    # print("protein = {}".format(len(protein)))
    # print("water = {}".format(len(water)))
    # print("ions = {}".format(len(ions)))
    # print("excipients = {}".format(len(excipients)))
    #
    # print("system = {}".format(len(system)))
    try:
        ls = expr.parseString(res_mask, parseAll=True)
    except Exception as e:
        raise ValueError("The passed restarint mask is not valid: {}".format(str(e)))

    atom_set = build_set(ls[0], dic_set)

    return atom_set
    
