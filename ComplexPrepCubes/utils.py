from openeye import oechem
from simtk.openmm import Vec3
from simtk import unit
from simtk.openmm import app
import itertools
from pkg_resources import resource_filename
from pdbfixer import PDBFixer
import parmed
from LigPrepCubes import ff_utils
from OpenMMCubes import utils
import logging

proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS',
                   'LEU', 'MET', 'PRO', 'THR', 'TYR',
                   'ARG', 'ASP', 'GLN', 'GLY', 'ILE',
                   'LYS', 'PHE', 'SER', 'TRP', 'VAL']

rnaResidues = ['A', 'G', 'C', 'U', 'I']
dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']


def split(mol):
    """
    This function splits the passed molecule in protein, ligand
    and water

    Parameters:
    ----------
    mol : oechem.OEMol
        The bio-molecular system to split

    Output:
    -------
    protein : oechem.OEMol
        The split protein
    ligand : oechem.OEMol
        The split ligand
    wat : oechel.OEMol
        The spit water which contenis the water and other
        excipient molecules
    
    """

    # Set empty molecule containers
    prot = oechem.OEMol()
    lig = oechem.OEMol()
    wat = oechem.OEMol()
    other = oechem.OEMol()

    # Define the Filter options before to split
    opt = oechem.OESplitMolComplexOptions()

    # The protein filter is set to avoid that multiple
    # chains are separated during the splitting
    pf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Protein)
    # The ligand filter is set to recognize just the ligand
    lf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Ligand)
    # The water filter is set to recognize just water molecules
    wf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Water)
    opt.SetProteinFilter(pf)
    opt.SetLigandFilter(lf)
    opt.SetWaterFilter(wf)
    
    # Splitting the system
    if not oechem.OESplitMolComplex(lig, prot, wat, other, mol, opt):
        oechem.OEThrow.Fatal('Unable to split the complex')
    
    # At this point prot contains the protein, lig contains the ligand,
    # wat contains the water and excipients contains the excipients

    return prot, lig, wat, other


def delete_shell(core_mol, del_mol, cut_off, in_out='in'):
    """
    This function deletes molecules present in the passed argument
    del_mol that are far (in_out=out) or close (in_out=in) than the
    selected cutoff distance (in A) from the passed molecules core_mol
       
    Parameters:
    -----------
    core_mol: OEMol molecule
        The core molecules 
    del_mol: OEMol molecule
        The molecules to be deleted if their distances from the core_mol
        molecules are greater or closer that the selected cutoff distance
    cut_off: python float number
        The threshold distance in A used to mark atom for deletion
    in_out: python string
        A flag used to select if delete molecules far or close than 
        the cutoff distance from the core_mol
    
    Return:
    -------
    to_del: copy of del_mol where atoms have been deleted 
    """

    if in_out not in ['in', 'out']:
        raise ValueError("The passed in_out parameter is not recognized: {}".format(in_out))

    # Copy the passed molecule to delete in
    to_del = oechem.OEMol(del_mol)

    # Create a OE bit vector mask for each atoms of the
    # molecule to delete
    bv = oechem.OEBitVector(to_del.GetMaxAtomIdx())
    bv.NegateBits()

    # Create the Nearest neighbours
    nn = oechem.OENearestNbrs(to_del, cut_off)
    for nbrs in nn.GetNbrs(core_mol):
        #bv.SetBitOff(nbrs.GetBgn().GetIdx())
        for atom in oechem.OEGetResidueAtoms(nbrs.GetBgn()):
            bv.SetBitOff(atom.GetIdx())

    # Invert selection mask
    if in_out == 'in':
        bv.NegateBits()

    pred = oechem.OEAtomIdxSelected(bv)
    for atom in to_del.GetAtoms(pred):
        to_del.DeleteAtom(atom)

    return to_del


def solvate(system, opt):
    # Load a fake system to initialize PDBfixer
    filename = resource_filename('pdbfixer', 'tests/data/test.pdb')
    fixer = PDBFixer(filename=filename)

    # Convert between OE and OpenMM topology
    omm_top, omm_pos = oemol_to_openmmTop(system)

    chain_names = []

    for chain in omm_top.chains():
        chain_names.append(chain.id)

    # Set the correct topology to the fake system
    fixer.topology = omm_top
    fixer.positions = omm_pos

    # Solvate the system
    fixer.addSolvent(padding=unit.Quantity(opt['solvent_padding'], unit.angstroms),
                     ionicStrength=unit.Quantity(opt['salt_concentration'], unit.millimolar))

    # The OpenMM topology produced by the solvation fixer has missing bond
    # orders and aromaticity. The following section is creating a new openmm
    # topology made of just water molecules and ions. The new topology is then
    # converted in an OEMol and added to the passed molecule to produce the
    # solvated system

    wat_ion_top = app.Topology()

    # Atom dictionary between the the fixer topology and the water_ion topology
    fixer_atom_to_wat_ion_atom = {}

    for chain in fixer.topology.chains():
        if chain.id not in chain_names:
            n_chain = wat_ion_top.addChain(chain.id)
            for res in chain.residues():
                n_res = wat_ion_top.addResidue(res.name, n_chain)
                for at in res.atoms():
                    n_at = wat_ion_top.addAtom(at.name, at.element, n_res)
                    fixer_atom_to_wat_ion_atom[at] = n_at

    for bond in fixer.topology.bonds():
        at0 = bond[0]
        at1 = bond[1]
        try:
            wat_ion_top.addBond(fixer_atom_to_wat_ion_atom[at0], fixer_atom_to_wat_ion_atom[at1], type=None, order=1)
        except:
            pass

    wat_ion_pos = fixer.positions[len(omm_pos):]

    oe_mol = openmmTop_to_oemol(wat_ion_top, wat_ion_pos)

    # Setting the box vectors
    omm_box_vectors = fixer.topology.getPeriodicBoxVectors()
    box_vectors = utils.PackageOEMol.encodePyObj(omm_box_vectors)
    oe_mol.SetData(oechem.OEGetTag('box_vectors'), box_vectors)

    oechem.OEAddMols(oe_mol, system)

    return oe_mol


def applyffProtein(protein, opt):

    # ofs = oechem.oemolostream("protein.oeb")
    # oechem.OEWriteMolecule(ofs, protein)

    # import sys
    # sys.exit(-1)

    topology, positions = oemol_to_openmmTop(protein)

    forcefield = app.ForceField(opt['protein_forcefield'])
    unmatched_residues = forcefield.getUnmatchedResidues(topology)

    if unmatched_residues:
        raise ValueError("The following Protein residues are not recognized "
                         "by the selected force field {}: {}".format(opt['protein_forcefield'], unmatched_residues))

    omm_system = forcefield.createSystem(topology, rigidWater=False)
    protein_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)

    return protein_structure


def applyffWater(water, opt):

    topology, positions = oemol_to_openmmTop(water)

    forcefield = app.ForceField(opt['solvent_forcefield'])
    unmatched_residues = forcefield.getUnmatchedResidues(topology)

    if unmatched_residues:
        raise ValueError("The following water molecules are not recognized "
                         "by the selected force field {}: {}".format(opt['solvent_forcefield'], unmatched_residues))

    omm_system = forcefield.createSystem(topology, rigidWater=False)
    water_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)

    return water_structure


def applyffExcipients(excipients, opt):

    topology, positions = oemol_to_openmmTop(excipients)

    forcefield = app.ForceField(opt['protein_forcefield'])
    unmatched_residues = forcefield.getUnmatchedResidues(topology)

    if unmatched_residues:
        raise ValueError("The following water molecules are not recognized "
                         "by the selected force field {}: {}".format(opt['protein_forcefield'], unmatched_residues))

    omm_system = forcefield.createSystem(topology, rigidWater=False)
    excipients_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)

    return excipients_structure


def applyffLigand(ligand, opt):

    # Check TLeap
    if opt['ligand_forcefield'] in ['GAFF', 'GAFF2']:
        ff_utils.ParamLigStructure(oechem.OEMol(), opt['ligand_forcefield']).checkTleap

    # Parametrize the Ligand
    pmd = ff_utils.ParamLigStructure(ligand, opt['ligand_forcefield'])
    molecule_structure = pmd.parameterize()
    molecule_structure.residues[0].name = "LIG"

    return molecule_structure


def oemol_to_openmmTop(mol):
    """
    This function converts an OEMol to an openmm topology
    The OEMol coordinates are assumed to be in Angstrom unit

    Parameters:
    -----------
    mol: OEMol molecule
        The molecule to convert 

    Return:
    -------
    topology : OpenMM Topology 
        The generated OpenMM topology
    positions : OpenMM Quantity
        The molecule atom positions associated wuth the 
        generated topology in Angstrom units
    """
    # OE Hierarchical molecule view
    hv = oechem.OEHierView(mol, oechem.OEAssumption_BondedResidue + oechem.OEAssumption_ResPerceived)

    # Create empty OpenMM Topology
    topology = app.Topology()
    # Dictionary used to map oe atoms to openmm atoms
    oe_atom_to_openmm_at = {}

    for chain in hv.GetChains():

        # Create empty OpenMM Chain
        openmm_chain = topology.addChain(chain.GetChainID())

        for frag in chain.GetFragments():

            for hres in frag.GetResidues():

                # Get OE residue
                oe_res = hres.GetOEResidue()
                # Create OpenMM residue
                openmm_res = topology.addResidue(oe_res.GetName(), openmm_chain)

                for oe_at in hres.GetAtoms():
                    # Select atom element based on the atomic number
                    element = app.element.Element.getByAtomicNumber(oe_at.GetAtomicNum())
                    # Add atom OpenMM atom to the topology
                    openmm_at = topology.addAtom(oe_at.GetName(), element, openmm_res)
                    openmm_at.index = oe_at.GetIdx()
                    # Add atom to the mapping dictionary
                    oe_atom_to_openmm_at[oe_at] = openmm_at

    # Create bonds preserving the bond ordering
    for bond in mol.GetBonds():
        aromatic = None

        # Set the bond aromaticity
        if bond.IsAromatic():
            aromatic = 'Aromatic'

        topology.addBond(oe_atom_to_openmm_at[bond.GetBgn()], oe_atom_to_openmm_at[bond.GetEnd()],
                         type=aromatic, order=bond.GetOrder())

    dic = mol.GetCoords()
    positions = [Vec3(v[0], v[1], v[2]) for k, v in dic.items()] * unit.angstrom

    return topology, positions


def openmmTop_to_oemol(topology, positions):
    """
    This function converts an OpenMM topology in an OEMol

    Parameters:
    -----------
    topology : OpenMM Topology 
        The OpenMM topology
    positions : OpenMM Quantity
        The molecule atom positions associated with the 
        topology


    Return:
    -------
    oe_mol : OEMol 
        The generated OEMol molecule
    """

    # Create an empty OEMol
    oe_mol = oechem.OEMol()

    # Mapping dictionary between openmm atoms and oe atoms
    openmm_atom_to_oe_atom = {}

    # Python set used to identify atoms that are not in protein residues
    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues)

    for chain in topology.chains():
        for res in chain.residues():
            # Create an OEResidue
            oe_res = oechem.OEResidue()
            # Set OEResidue name
            oe_res.SetName(res.name)
            # If the atom is not a protein atom then set its heteroatom
            # flag to True
            if res.name not in keep:
                oe_res.SetFragmentNumber(chain.index + 1)
                oe_res.SetHetAtom(True)
            # Set OEResidue Chain ID
            oe_res.SetChainID(chain.id)
            # res_idx = int(res.id) - chain.index * len(chain._residues)
            # Set OEResidue number
            oe_res.SetResidueNumber(int(res.id))

            for openmm_at in res.atoms():
                # Create an OEAtom  based on the atomic number
                oe_atom = oe_mol.NewAtom(openmm_at.element._atomic_number)
                # Set atom name
                oe_atom.SetName(openmm_at.name)
                # Set Atom index
                oe_res.SetSerialNumber(openmm_at.index + 1)
                # Commit the changes
                oechem.OEAtomSetResidue(oe_atom, oe_res)
                # Update the dictionary OpenMM to OE
                openmm_atom_to_oe_atom[openmm_at] = oe_atom

    # Create the bonds
    for bond in topology.bonds():
        at0 = bond[0]
        at1 = bond[1]
        # Read in the bond order from the OpenMM topology
        bond_order = bond.order

        # If bond order info are not present set the bond order to one
        if not bond_order:
            logging.info("WARNING: Bond order info missing between atom indexes: {}-{}".format(at0.index, at1.index))
            bond_order = 1

        # OE atoms
        oe_atom0 = openmm_atom_to_oe_atom[at0]
        oe_atom1 = openmm_atom_to_oe_atom[at1]

        # Set OE atom aromaticity
        if bond.type:
            oe_atom0.SetAromatic(True)
            oe_atom1.SetAromatic(True)

        # Create the bond
        oe_bond = oe_mol.NewBond(oe_atom0, oe_atom1, bond_order)

        if bond.type:
            oe_bond.SetAromatic(True)

    # Set the OEMol positions
    pos = positions.in_units_of(unit.angstrom) / unit.angstrom
    pos = list(itertools.chain.from_iterable(pos))
    oe_mol.SetCoords(pos)

    return oe_mol


def order_check(mol, fname):
    import logging
    logger = logging.getLogger('Testing')
    hdlr = logging.FileHandler(fname)
    formatter = logging.Formatter('%(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)

    # hv = oechem.OEHierView(mol, oechem.OEAssumption_BondedResidue +
    #                        oechem.OEAssumption_ResPerceived)
    hv = oechem.OEHierView(mol)


    for chain in hv.GetChains():
        logger.info('{}'.format(chain.GetChainID()))
        for frag in chain.GetFragments():
            for hres in frag.GetResidues():
                logger.info('\t{} {}'.format(hres.GetOEResidue().GetName(), hres.GetOEResidue().GetResidueNumber()))
                for oe_at in hres.GetAtoms():
                    logger.info('\t\t{}'.format(oe_at.GetName()))

    return