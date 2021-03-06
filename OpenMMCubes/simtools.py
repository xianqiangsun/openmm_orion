import sys, mdtraj, tarfile, os
import numpy as np
from sys import stdout
from openeye import oechem
from simtk import unit, openmm
from simtk.openmm import app
import pyparsing as pyp
from floe.api.orion import in_orion,  upload_file
from oeommtools import utils as oeommutils

try:
    import cPickle as pickle
except ImportError:
    import pickle


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
            delta = box_v/2 - cog
            # New Coordinates
            new_coords = coords + delta
            structure.coordinates = new_coords
            positions = structure.positions

    # OpenMM system
    system = structure.createSystem(nonbondedMethod=eval("app.%s" % opt['nonbondedMethod']),
                                    nonbondedCutoff=opt['nonbondedCutoff']*unit.angstroms,
                                    constraints=eval("app.%s" % opt['constraints']),
                                    removeCMMotion=False)
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
        try:
            platform = openmm.Platform.getPlatformByName(opt['platform'])
        except Exception as e:
            oechem.OEThrow.Fatal('The selected platform is not supported: {}'.format(str(e)))

        try:
            # Set platform CUDA or OpenCL precision
            properties = {'Precision': opt['cuda_opencl_precision']}

            simulation = app.Simulation(topology, system, integrator,
                                        platform=platform,
                                        platformProperties=properties)
        except Exception:
            oechem.OEThrow.Fatal('It was not possible to set the {} precision for the {} platform'
                                 .format(opt['cuda_opencl_precision'], opt['platform']))

    # Set starting positions and velocities
    simulation.context.setPositions(positions)

    # Set Box dimensions
    simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

    # If the velocities are not present in the Parmed structure
    # new velocity vectors are generated otherwise the system is
    # restarted from the previous State
    if opt['SimType'] in ['nvt', 'npt']:

        if opt['trajectory_interval']:
            structure.save(opt['outfname']+'.pdb', overwrite=True)
            # GAC ADDED - TESTING
            # Preserve original pdb file residue numbers
            pdbfname_test = opt['outfname'] + '_ordering_test' + '.pdb'
            ofs = oechem.oemolostream(pdbfname_test)
            flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms
            ofs.SetFlavor(oechem.OEFormat_PDB, flavor)

            new_temp_mol = oeommutils.openmmTop_to_oemol(structure.topology, structure.positions, verbose=False)
            new_pos = new_temp_mol.GetCoords()
            opt['molecule'].SetCoords(new_pos)
            oechem.OEWriteConstMolecule(ofs, opt['molecule'])

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
        for k, v in uname()._asdict().items():
            print(k, ':', v, file=printfile)
        # Platform properties
        for prop in mmplat.getPropertyNames():
            val = mmplat.getPropertyValue(simulation.context, prop)
            print(prop, ':', val, file=printfile)
            
    print('OpenMM({}) simulation generated for {} platform'.format(mmver, mmplat.getName()), file=printfile)

    if opt['SimType'] in ['nvt', 'npt']:

        opt['Logger'].info('Running {time} ps = {steps} steps of {SimType} at {temperature} K'.format(**opt))
        
        # Start Simulation
        simulation.step(opt['steps'])

        state = simulation.context.getState(getPositions=True, getVelocities=True,
                                            getEnergy=True, enforcePeriodicBox=True)
        
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

        # If required uploading files to Orion
        _file_processing(**opt)

    # Update the OEMol complex positions to match the new
    # Parmed structure after the simulation
    new_temp_mol = oeommutils.openmmTop_to_oemol(structure.topology, structure.positions, verbose=False)
    new_pos = new_temp_mol.GetCoords()
    opt['molecule'].SetCoords(new_pos)

    return


def _file_processing(**opt):
    """
    This supporting function compresses the produced trajectory
    and supporting files in a .tar file (if required ) and eventually
    uploaded them to Orion. If not .tar file is selected then all the
    generated files are eventually uploaded in Orion

    Parameters
    ----------
    opt: python dictionary
        A dictionary containing all the MD setting info
    """

    # Set the trajectory file name
    if opt['trajectory_filetype'] == 'NetCDF':
        trj_fn = opt['outfname'] +'.nc'
    elif opt['trajectory_filetype'] == 'DCD':
        trj_fn = opt['outfname'] +'.dcd'
    elif opt['trajectory_filetype'] == 'HDF5':
        trj_fn = opt['outfname'] + '.hdf5'
    else:
        oechem.OEThrow.Fatal("The selected trajectory filetype is not supported: {}"
                             .format(opt['trajectory_filetype']))
    # Set .pdb file names
    pdb_fn = opt['outfname'] + '.pdb'
    pdb_order_fn = opt['outfname'] + '_ordering_test' + '.pdb'
    log_fn = opt['outfname'] + '.log'

    # List all the file names
    fnames = [trj_fn, pdb_fn, pdb_order_fn, log_fn]

    ex_files = []

    # Check which file names are actually produced files
    for fn in fnames:
        if os.path.isfile(fn):
            ex_files.append(fn)

    # Tar the outputted files if required
    if opt['tarxz']:

        tarname = opt['outfname'] + '.tar.xz'

        opt['Logger'].info('Creating tarxz file: {}'.format(tarname))

        tar = tarfile.open(tarname, "w:xz")

        for name in ex_files:
            opt['Logger'].info('Adding {} to {}'.format(name, tarname))
            tar.add(name)
        tar.close()

        if in_orion():
            upload_file(tarname, tarname, tags=['TRJ_INFO'])

        # Clean up files that have been added to tar.
        for tmp in ex_files:
            try:
                os.remove(tmp)
            except:
                pass
    else:  # If not .tar file is required the files are eventually uploaded in Orion
        if in_orion():
            for fn in ex_files:
                upload_file(fn, fn, tags=['TRJ_INFO'])

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

    reporters = []

    if opt['reporter_interval']:
        state_reporter = app.StateDataReporter(outfname+'.log', separator="\t",
                                               reportInterval=opt['reporter_interval'],
                                               step=True,
                                               potentialEnergy=True, totalEnergy=True,
                                               volume=True, temperature=True)

        reporters.append(state_reporter)

        progress_reporter = app.StateDataReporter(stdout, separator="\t",
                                                  reportInterval=opt['reporter_interval'],
                                                  step=True, totalSteps=totalSteps,
                                                  time=True, speed=True, progress=True,
                                                  elapsedTime=True, remainingTime=True)

        reporters.append(progress_reporter)

    if opt['trajectory_interval']:

        # Trajectory file format selection
        if opt['trajectory_filetype'] == 'NetCDF':
            traj_reporter = mdtraj.reporters.NetCDFReporter(outfname+'.nc', opt['trajectory_interval'])
        elif opt['trajectory_filetype'] == 'DCD':
            traj_reporter = app.DCDReporter(outfname+'.dcd', opt['trajectory_interval'])
        elif opt['trajectory_filetype'] == 'HDF5':
            mdtraj.reporters.HDF5Reporter(outfname+'.hdf5', opt['trajectory_interval'])
        else:
            oechem.OEThrow.Fatal("The selected trajectory file format is not supported: {}"
                                 .format(opt['trajectory_filetype']))

        reporters.append(traj_reporter)

    return reporters


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

    def split(system):
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
        # cofactor_set = set()
        # system_set = set()

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
