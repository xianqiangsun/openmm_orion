from openeye import oechem
from simtk.openmm import Vec3
from simtk import unit
from simtk.openmm import app
import itertools


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
    # chains are separeted during the splitting
    pf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Protein)
    opt.SetProteinFilter(pf)

    
    # Splitting the system
    oechem.OESplitMolComplex(lig, prot, wat, other, mol, opt)
    
    # At this point prot contains the protein, lig contains the ligand and
    # wat contains the water and excipients

    return prot, wat, other


def delete_shell(core_mol, del_mol, cut_off, in_out='in'):
    """
    This function keeps the molecules in del_mol that are 
    inside or outside (in_out) a shell built around the 
    core_mol. The cutoff distance (in A) select the shell 
    threshold distance between the core_mol and the 
    del_mol atoms.
      
    Parameters:
    -----------
    core_mol: OEMol molecule
        The core molecule 
    del_mol: OEMol molecule
        The molecule to delete atoms inside or outside the shell
    cut_off: python float number
        The threshold distance in A used to mark atom for deletion
    in_out: python string
        A flag used to select which atoms to delete. 'in' and 'out' 
        will respectively delete all the atoms that are inside or 
        outside the shell.
    
    Return:
    -------
    to_del: the copy of del_mol where atoms have been deleted 
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


def oemol_to_openmmTop(mol):

    hv = oechem.OEHierView(mol)

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

    # Create bonds
    for bond in mol.GetBonds():
        topology.addBond(oe_atom_to_openmm_at[bond.GetBgn()], oe_atom_to_openmm_at[bond.GetEnd()])

    dic = mol.GetCoords()
    positions = [Vec3(v[0], v[1], v[2]) for k, v in dic.items()] * unit.angstrom

    return topology, positions


def openmmTop_to_oemol(topology, positions):

    oe_mol = oechem.OEMol()
    openmm_atom_to_oe_atom={}

    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues)

    for chain in topology.chains():
        for res in chain.residues():
            oe_res = oechem.OEResidue()
            oe_res.SetName(res.name)
            if res.name not in keep:
                oe_res.SetHetAtom(True)
            oe_res.SetChainID(chain.id)
            #res_idx = int(res.id) - chain.index * len(chain._residues)
            oe_res.SetResidueNumber(int(res.id))
            for openmm_at in res.atoms():
                oe_atom = oe_mol.NewAtom(openmm_at.element._atomic_number)
                oe_atom.SetName(openmm_at.name)
                oe_res.SetSerialNumber(openmm_at.index+1)
                oechem.OEAtomSetResidue(oe_atom, oe_res)

                openmm_atom_to_oe_atom[openmm_at] = oe_atom

    for bond in topology.bonds():
        at0 = bond[0]
        at1 = bond[1]
        oe_mol.NewBond(openmm_atom_to_oe_atom[at0], openmm_atom_to_oe_atom[at1], 1)

    #oe_coords = oechem.OEFloatArray(oe_mol.GetMaxAtomIdx() * 3)

    pos = positions.in_units_of(unit.angstrom)/unit.angstrom

    pos = list(itertools.chain.from_iterable(pos))

    oe_mol.SetCoords(pos)

    return oe_mol
