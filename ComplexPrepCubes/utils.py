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
    oechem.OESplitMolComplex(lig, prot, wat, other, mol, opt)
    
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

    # Create bonds preserving the bond ordering assessed by the OE Toolkit
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
        The molecule atom positions associated wuth the 
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
            raise RuntimeWarning("Bond order info missing between atom indexes: {}-{}".format(at0.index, at1.index))
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





