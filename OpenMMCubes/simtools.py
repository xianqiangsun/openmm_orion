import io, os, sys, base64, parmed, mdtraj, pdbfixer
import numpy as np
from sys import stdout
from openeye import oechem
from simtk import unit, openmm
from simtk.openmm import app
import OpenMMCubes.utils as utils
import pyparsing as pyp
from ComplexPrepCubes import utils as prep_utils

try:
    import cPickle as pickle
except ImportError:
    import pickle

    
# def genProteinStructure(proteinpdb, **opt):
#     """
#     Starting from OpenMM PDBFile, generates the OpenMM System for the protein and
#     then create the parametrized ParmEd Structure of the protein.
#
#     Parameters
#     ----------
#     proteinpdb : openmm.app.PDBFile object,
#         Loaded PDBFile object of the protein.
#     protein_forcefield : (opt), xml file, default='amber99sbildn.xml'
#         Forcefield parameters for protein
#     solvent_forcefield : opt), xml file, default='tip3p.xml'
#         Forcefield parameters for solvent
#
#     Returns
#     -------
#     solv_structure : parmed.structure.Structure
#         The parameterized Structure of the protein with solvent molecules. (No ligand).
#     """
#     #Generate protein Structure object
#     forcefield = app.ForceField(opt['protein_forcefield'], opt['solvent_forcefield'])
#     protein_system = forcefield.createSystem( proteinpdb.topology )
#     protein_structure = parmed.openmm.load_topology(proteinpdb.topology,
#                                                     protein_system,
#                                                     xyz=proteinpdb.positions)
#     return protein_structure
#
#
# def combinePositions(proteinPositions, molPositions):
#     """
#     Loops through the positions from the ParmEd structures of the protein and ligand,
#     divides by unit.angstroms which will ensure both positions arrays are in the same units.
#
#     Parameters
#     ----------
#     proteinPositions : list of 3-element Quantity tuples.
#         Positions list taken directly from the protein Structure.
#     molPositions : list of 3-element Quantity tuples.
#         Positions list taken directly from the molecule Structure.
#
#     Returns
#     -------
#     positions : list of 3-element Quantity tuples.
#         ex. unit.Quantity(positions, positions_unit)
#         Combined positions of the protein and molecule Structures.
#     """
#     positions_unit = unit.angstroms
#     positions0_dimensionless = np.array(proteinPositions / positions_unit)
#     positions1_dimensionless = np.array(molPositions / positions_unit)
#     coordinates = np.vstack(
#         (positions0_dimensionless, positions1_dimensionless))
#     natoms = len(coordinates)
#     positions = np.zeros([natoms, 3], np.float32)
#     for index in range(natoms):
#             (x, y, z) = coordinates[index]
#             positions[index, 0] = x
#             positions[index, 1] = y
#             positions[index, 2] = z
#     positions = unit.Quantity(positions, positions_unit)
#     return positions
#
#
# def mergeStructure(proteinStructure, molStructure):
#     """
#     Combines the parametrized ParmEd structures of the protein and ligand to
#     create the Structure for the protein:ligand complex, while retaining the SMIRFF
#     parameters on the ligand. Preserves positions and box vectors.
#     (Not as easily achieved using native OpenMM tools).
#
#     Parameters
#     ----------
#     proteinStructure : parmed.structure.Structure
#         The parametrized structure of the protein.
#     moleculeStructure : parmed.structure.Structure
#         The parametrized structure of the ligand.
#
#     Returns
#     -------
#     structure : parmed.structure.Structure
#         The parametrized structure of the protein:ligand complex.
#     """
#     structure = proteinStructure + molStructure
#     positions = combinePositions(proteinStructure.positions, molStructure.positions)
#     # Concatenate positions arrays (ensures same units)
#     structure.positions = positions
#     # Restore original box vectors
#     structure.box = proteinStructure.box
#     return structure
#
#
# def solvateComplexStructure(structure, **opt):
#     """
#     Uses PDBFixer to add missing atoms, assign protonation states, and solvate
#     the parametrized ParmEd Structure of the protein:ligand complex. Returns the
#     parametrized Structure without the ligand present.
#     We initially combine the ligand with the protein to prevent PDBFixer
#     from adding solvent into the binding pocket.
#
#     Parameters
#     ----------
#     structure : parmed.structure.Structure
#         The parametrized structure of the protein:ligand complex
#     pH : (opt), float, default=7.4
#         Solvent pH used to select appropriate protein protonation state.
#     solvent_padding : (opt), float, default=10
#         Padding around protein for solvent box (unit.angstroms)
#     salt_concentration : (opt), float, default=10
#         Salt concentration (millimolar)
#     protein_forcefield : (opt), xml file, default='amber99sbildn9sbildn.xml'
#         Forcefield parameters for protein
#     solvent_forcefield : opt), xml file, default='tip3p.xml'
#         Forcefield parameters for solvent
#
#     Returns
#     -------
#     solv_structure : parmed.structure.Structure
#         The parameterized Structure of the protein with solvent molecules. (No ligand).
#     """
#     if opt['Logger'] is None:
#         printfile = sys.stdout
#     else:
#         printfile = opt['Logger'].file
#
#     tmpfile = opt['outfname']+'-pl.tmp'
#     structure.save(tmpfile,format='pdb')
#
#     seqres = False
#     with open(tmpfile, 'r') as infile:
#         for line in infile:
#             if 'SEQRES' in line:
#                 seqres = True
#                 break
#     if not seqres:
#         print('Did not find SEQRES in PDB. PDBFixer will not find missing Residues.', file=printfile)
#
#     # Solvate with PDBFixer
#     print('PDBFixer settings for {}'.format(opt['outfname']), file=printfile)
#     print('\tpH = {}'.format(opt['pH']), file=printfile)
#     print('\tpadding = {}'.format(unit.Quantity(opt['solvent_padding'], unit.angstroms)), file=printfile)
#     print('\tsalt conc. = {}'.format(unit.Quantity(opt['salt_concentration'], unit.millimolar)), file=printfile)
#     fixer = pdbfixer.PDBFixer(tmpfile)
#     fixer.findMissingResidues()
#     fixer.findNonstandardResidues()
#     fixer.findMissingAtoms()
#
#     if fixer.missingAtoms:
#         print('Found missing Atoms:', fixer.missingAtoms, file=printfile)
#     if fixer.missingTerminals:
#         print('Found missing Terminals:', fixer.missingTerminals, file=printfile)
#     if fixer.nonstandardResidues:
#         print('Found nonstandard Residues:', fixer.nonstandardResidues, file=printfile)
#
#     fixer.replaceNonstandardResidues()
#     #fixer.removeHeterogens(False)
#     fixer.addMissingAtoms()
#     fixer.addMissingHydrogens(opt['pH'])
#     fixer.addSolvent(padding=unit.Quantity(opt['solvent_padding'], unit.angstroms),
#                 ionicStrength=unit.Quantity(opt['salt_concentration'], unit.millimolar))
#
#     # Load PDBFixer object back to Structure
#     tmp = parmed.openmm.load_topology(fixer.topology, xyz=fixer.positions)
#
#     # Remove ligand from protein Structure by AmberMask selection
#     tmp.strip(":LIG")
#     tmp.save(opt['outfname']+'-nomol.tmp',format='pdb')
#     # Reload PDBFile
#     nomol = app.PDBFile(opt['outfname']+'-nomol.tmp')
#     forcefield = app.ForceField(opt['protein_forcefield'], opt['solvent_forcefield'])
#     nomol_system = forcefield.createSystem(nomol.topology, rigidWater=False)
#     # Regenerate parameterized solvated protein structure
#     solv_structure = parmed.openmm.load_topology(nomol.topology,
#                                                 nomol_system,
#                                                 xyz=nomol.positions)
#     # Restore box vectors
#     solv_structure.box = tmp.box
#
#     tmpfiles = [ opt['outfname']+'-pl.tmp', opt['outfname']+'-nomol.tmp' ]
#     utils.cleanup(tmpfiles)
#
#     return solv_structure


def simulation(mdData, **opt):
    """
    This supporting function performs: OpenMM Minimization, NVT and NPT
    Molecular Dynamics (MD) simulations

    Parameters
    ----------
    mdData : MDData data object
        The object which recovers the relevant Parmed structure data
        to perform MD
    opt: python dictionary
        A dictionary containing all the MD setting info
    """
    
    if opt['Logger'] is None:
        printfile = sys.stdout
    else:
        printfile = opt['Logger'].file

    # MD data extracted from Parmed
    structure = mdData.structure
    topology = mdData.topology
    positions = mdData.positions
    velocities = mdData.velocities
    box = mdData.box

    # Time step in ps
    stepLen = 0.002 * unit.picoseconds

    if opt['SimType'] in ['nvt', 'npt']:
        # Centering the system to the OpenMM Unit Cell
        if opt['center']:
            opt['Logger'].info("Centering is On")
            # Numpy array in A
            coords = structure.coordinates
            # System Center of Geometry
            cog = np.mean(coords, axis=0)
            # System box vectors
            box_v = structure.box_vectors.in_units_of(unit.angstrom)/unit.angstrom
            box_v = np.array([box_v[0][0], box_v[1][1], box_v[2][2]])
            # Translation vector
            delta = box_v/2 -cog
            # New Coordinates
            new_coords = coords + delta
            structure.coordinates = new_coords
            positions = structure.positions

    # OpenMM system
    system = structure.createSystem(nonbondedMethod=eval("app.%s" % opt['nonbondedMethod']),
                                    nonbondedCutoff=opt['nonbondedCutoff']*unit.angstroms,
                                    constraints=eval("app.%s" % opt['constraints']))
    # OpenMM Integrator
    integrator = openmm.LangevinIntegrator(opt['temperature']*unit.kelvin, 1/unit.picoseconds, stepLen)
    
    if opt['SimType'] == 'npt':
        # Add Force Barostat to the system
        system.addForce(openmm.MonteCarloBarostat(opt['pressure']*unit.atmospheres, opt['temperature']*unit.kelvin, 25))

    # Apply restraints
    if opt['restraints']:
        opt['Logger'].info("RESTRAINT mask applied to: {}"
                           "\tRestraint weight: {}".format(opt['restraints'],
                                                           opt['restraintWt'] *
                                                           unit.kilocalories_per_mole/unit.angstroms**2))
        # Select atom to restraint
        res_atom_set = restraints(opt['molecule'], mask=opt['restraints'])
        opt['Logger'].info("Number of restraint atoms: {}".format(len(res_atom_set)))
        # define the custom force to restrain atoms to their starting positions
        force_restr = openmm.CustomExternalForce('k_restr*periodicdistance(x, y, z, x0, y0, z0)^2')
        # Add the restraint weight as a global parameter in kcal/mol/A^2
        force_restr.addGlobalParameter("k_restr", opt['restraintWt']*unit.kilocalories_per_mole/unit.angstroms**2)
        # Define the target xyz coords for the restraint as per-atom (per-particle) parameters
        force_restr.addPerParticleParameter("x0")
        force_restr.addPerParticleParameter("y0")
        force_restr.addPerParticleParameter("z0")

        for idx in range(0, len(positions)):
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
    
    # Set starting positions and velocities
    simulation.context.setPositions(positions)

    # Set Box dimensions
    simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

    # If the velocities are not present in the Parmed structure
    # new velocity vectors are generated otherwise the system is
    # restarted from the previous State
    if opt['SimType'] in ['nvt', 'npt']:
        if velocities is not None:
            opt['Logger'].info('RESTARTING simulation from a previous State')
            simulation.context.setVelocities(velocities)
        else:
            # Set the velocities drawing from the Boltzmann distribution at the selected temperature
            opt['Logger'].info('GENERATING a new starting State')
            simulation.context.setVelocitiesToTemperature(opt['temperature']*unit.kelvin)

        # Convert simulation time in steps
        opt['steps'] = int(round(opt['time']/(stepLen.in_units_of(unit.picoseconds)/unit.picoseconds)))
        
        # Set Reporters
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

        opt['Logger'].info('Running {time} ps = {steps} steps of {SimType} MD at {temperature} K'.format( **opt))
        
        # Start Simulation
        simulation.step(opt['steps'])

        if opt['convert']:
            opt['Logger'].info('Converting trajectories to: {trajectory_filetype}'.format(**opt))
            mdTrajConvert(simulation, outfname=opt['outfname'],
                               trajectory_selection=opt['trajectory_selection'],
                               trajectory_filetype=opt['trajectory_filetype'])

        state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True)
        
    elif opt['SimType'] == 'min':
        # Start Simulation
        opt['Logger'].info('Minimization steps: {steps}'.format(**opt))

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


def restraints(system, mask=''):
    """
    This function selects the atom indexes to apply harmonic restraints

    Parameters
    ----------
    system : OEMol with the attached Parmed structure
        The system to apply restraints to

    mask : python string
        A string used to select the atom to restraints. A Backus–Naur Form grammar
        (https://en.wikipedia.org/wiki/Backus–Naur_form) is defined by using pyparsing.
        The defined tokens are: "ligand", "protein", "ca_protein" ,"water", "ions",
        "cofactors" and "resid chain1:res_idx1 chain2:res_idx2 ... res_idxn"

        Operator tokens are:
        "not" = invert selection
        "or" = add selections
        "and" = intersect selections
        "diff" = logic difference between selection
        "noh" = remove hydrogens from the selection

    Returns
    -------
    atom_set : python set
        the select atom indexes to apply harmonic restraints

    Notes
    -----
        Example of selection string:
        mask = "ligand or protein"
        mask = "not water or not ions"
        mask = "ligand or protein or cofactors"
        mask = "noh protein"
        mask = "resid A:17 B:12 17 18"
        mask = "protein diff resid A:1"
    """

    def split(mol):
        """
        This function splits the passed molecule in components and tracks the
        mapping between the original molecule and the split components. The
        mapping is created as separated atom component index sets.

        Parameters:
        -----------
        mol: OEMol
            The molecule to split in components. The components are:
                the protein atoms,
                the protein carbon alpha atoms
                the water atoms,
                the ion atoms,
                the cofactor atoms
        Returns:
        --------
        dic_set: python dictionary
            The dictionary has token words as keys and for value the related
            atom set. The token keywords are:
                protein,
                ca_protein,
                ligand,
                water,
                ions,
                cofactors,
                system
        """

        # Define Empty sets
        lig_set = set()
        prot_set = set()
        ca_prot_set = set()
        wat_set = set()
        excp_set = set()
        ion_set = set()
        cofactor_set = set()
        system_set = set()

        # Atom Bond Set vector used to contains the whole system
        frags = oechem.OEAtomBondSetVector()

        # Define Options for the Filter
        opt = oechem.OESplitMolComplexOptions()

        # The protein filter is set to avoid that multiple
        # chains are separated during the splitting
        pf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Protein)

        # The ligand filter is set to recognize just the ligand
        lf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Ligand)

        # The water filter is set to recognize just water molecules
        wf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Water)

        # Set options based on the defined filters
        opt.SetProteinFilter(pf)
        opt.SetLigandFilter(lf)
        opt.SetWaterFilter(wf)

        # Define the system fragments
        if not oechem.OEGetMolComplexFragments(frags, system, opt):
            oechem.OEThrow.Fatal('Unable to generate the system fragments')

        # Set empty OEMol containers
        prot = oechem.OEMol()
        lig = oechem.OEMol()
        wat = oechem.OEMol()
        excp = oechem.OEMol()
        import logging
        # Split the protein from the system
        atommap = oechem.OEAtomArray(system.GetMaxAtomIdx())
        if not oechem.OECombineMolComplexFragments(prot, frags, opt, opt.GetProteinFilter(), atommap):
            oechem.OEThrow.Fatal('Unable to split the Protein')
        # Populate the protein set and the protein carbon alpha set
        pred = oechem.OEIsAlphaCarbon()
        for sys_at in system.GetAtoms():
            sys_idx = sys_at.GetIdx()
            at_idx = atommap[sys_idx]
            if at_idx:
                prot_set.add(sys_idx)
                at = system.GetAtom(oechem.OEHasAtomIdx(sys_idx))
                if pred(at):
                    ca_prot_set.add(sys_idx)
                # print(sys_idx, '->', at_idx)

        # Split the ligand from the system
        atommap = oechem.OEAtomArray(system.GetMaxAtomIdx())
        if not oechem.OECombineMolComplexFragments(lig, frags, opt, opt.GetLigandFilter(), atommap):
            oechem.OEThrow.Fatal('Unable to split the Ligand')
        # Populate the ligand set
        for sys_at in system.GetAtoms():
            sys_idx = sys_at.GetIdx()
            at_idx = atommap[sys_idx]
            if at_idx:
                lig_set.add(sys_idx)
                # print(sys_idx, '->', at_idx)

        # Split the water from the system
        atommap = oechem.OEAtomArray(system.GetMaxAtomIdx())
        if not oechem.OECombineMolComplexFragments(wat, frags, opt, opt.GetWaterFilter(), atommap):
            oechem.OEThrow.Fatal('Unable to split the Water')
        # Populate the water set
        for sys_at in system.GetAtoms():
            sys_idx = sys_at.GetIdx()
            at_idx = atommap[sys_idx]
            if at_idx:
                wat_set.add(sys_idx)
                # print(sys_idx, '->', at_idx)

        # Split the excipients from the system
        atommap = oechem.OEAtomArray(system.GetMaxAtomIdx())
        if not oechem.OECombineMolComplexFragments(excp, frags, opt, opt.GetOtherFilter(), atommap):
             oechem.OEThrow.Fatal('Unable to split the Excipients')
        # Populate the excipient set
        for sys_at in system.GetAtoms():
            sys_idx = sys_at.GetIdx()
            at_idx = atommap[sys_idx]
            if at_idx:
                excp_set.add(sys_idx)
                # print(sys_idx, '->', at_idx)

        # Create the ions set
        for exc_idx in excp_set:
            atom = system.GetAtom(oechem.OEHasAtomIdx(exc_idx))
            if atom.GetDegree() == 0:
                ion_set.add(exc_idx)

        # Create the cofactor set
        cofactor_set = excp_set - ion_set

        # Create the system set
        system_set = prot_set | lig_set | excp_set | wat_set

        if len(system_set) != system.NumAtoms():
            oechem.OEThrow.Fatal("The total system atom number {} is different "
                                 "from its set representation {}".format(system.NumAtoms(), system_set))

        # The dictionary is used to link the token keywords to the created molecule sets
        dic_set = {'ligand': lig_set, 'protein': prot_set, 'ca_protein': ca_prot_set,
                   'water': wat_set,  'ions': ion_set,     'cofactors': cofactor_set, 'system': system_set}

        return dic_set

    def build_set(ls, dsets):
        """
        This function select the atom indexes to apply restraints

        Parameters:
        -----------
        ls: python list
            the parsed list with tokens and operand tokes for the selection
        dsets: python dictionary
             the dictionary containing the sets for the selection

        Return:
        -------
        atom_set: python set
            the set containing the atom index to restrain
        """

        def noh(ls, dsets):
            """
            This function remove hydrogens from the selection
            """
            data_set = build_set(ls[1], dsets)

            noh_set = set()
            pred = oechem.OEIsHydrogen()

            for idx in data_set:
                atom = system.GetAtom(oechem.OEHasAtomIdx(idx))
                if not pred(atom):
                    noh_set.add(idx)

            return noh_set

        def residues(ls):
            """
            This function select residues based on the residue numbers. An example of
            selection can be:
            mask = 'resid A:16 17 19 B:1'
            """
            # List residue atom index to be restrained
            res_atom_set = set()

            # Dictionary of lists with the chain residues selected to be restrained
            # e.g. {chainA:[res1, res15], chainB:[res19, res17]}
            chain_dic = {'': []}

            # Fill out the chain dictionary
            i = 0
            while i < len(ls):
                if ls[i].isdigit():
                    chain_dic[''].append(int(ls[i]))
                    i += 1
                else:
                    try:
                        chain_dic[ls[i]].append(int(ls[i + 2]))
                    except:
                        chain_dic[ls[i]] = []
                        chain_dic[ls[i]].append(int(ls[i + 2]))
                    i += 3

            # Loop over the molecular system to select the atom indexes to be restrained
            hv = oechem.OEHierView(system, oechem.OEAssumption_BondedResidue + oechem.OEAssumption_ResPerceived)
            for chain in hv.GetChains():
                chain_id = chain.GetChainID()
                if chain_id not in chain_dic:
                    continue
                for frag in chain.GetFragments():
                    for hres in frag.GetResidues():
                        res_num = hres.GetOEResidue().GetResidueNumber()
                        if res_num not in chain_dic[chain_id]:
                            continue
                        for oe_at in hres.GetAtoms():
                            res_atom_set.add(oe_at.GetIdx())

            return res_atom_set

        # Terminal Literal return the related set
        if isinstance(ls, str):
            return dsets[ls]
        # Not or Noh
        if len(ls) == 2:
            if ls[0] == 'noh':  # Noh case
                return noh(ls, dsets)
            elif ls[0] == 'not':  # Not case
                return dsets['system'] - build_set(ls[1], dsets)
            else:  # Resid case with one index
                return residues(ls[1])

        if len(ls) == 3:
            if ls[1] == 'or':  # Or Case (set union)
                return build_set(ls[0], dsets) | build_set(ls[2], dsets)
            elif ls[1] == 'and':  # And Case (set intersection)
                return build_set(ls[0], dsets) & build_set(ls[2], dsets)
            elif ls[1] == 'diff':  # Diff case (set difference)
                return build_set(ls[0], dsets) - build_set(ls[2], dsets)
            else:
                return residues(ls[1:])  # Resid case with one or two indexes
        else:
            if ls[0] == 'resid':
                return residues(ls[1:])  # Resid case with multiple indexes
            else:
                raise ValueError("The passed list have too many tokens: {}".format(ls))

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

    # Restraint  function body

    # Residue number selection
    id = pyp.Optional(pyp.Word(pyp.alphanums) + pyp.Literal(':')) + pyp.Word(pyp.nums)
    resid = pyp.Group(pyp.Literal("resid") + pyp.OneOrMore(id))

    # Define the tokens for the BNF grammar
    operand = pyp.Literal("protein") | pyp.Literal("ca_protein") | \
              pyp.Literal("ligand") | pyp.Literal("water") | \
              pyp.Literal("ions") | pyp.Literal("cofactors") | resid

    # BNF Grammar definition with parseAction makeLRlike
    expr = pyp.operatorPrecedence(operand,
                                    [
                                        (None, 2, pyp.opAssoc.LEFT, makeLRlike(None)),
                                        (pyp.Literal("not"), 1, pyp.opAssoc.RIGHT, makeLRlike(1)),
                                        (pyp.Literal("noh"), 1, pyp.opAssoc.RIGHT, makeLRlike(1)),
                                        (pyp.Literal("and"), 2, pyp.opAssoc.LEFT, makeLRlike(2)),
                                        (pyp.Literal("or"), 2, pyp.opAssoc.LEFT, makeLRlike(2)),
                                        (pyp.Literal("diff"), 2, pyp.opAssoc.LEFT, makeLRlike(2))
                                    ])
    # Parse the input string
    try:
        ls = expr.parseString(mask, parseAll=True)
    except Exception as e:
        raise ValueError("The passed restraint mask is not valid: {}".format(str(e)))

    # Split the system
    dic_sets = split(system)

    # Select the atom index to restraints
    atom_set = build_set(ls[0], dic_sets)

    return atom_set
