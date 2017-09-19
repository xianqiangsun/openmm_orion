from openeye import oechem
from openeye import oequacpac
from simtk.openmm import app
import parmed
from LigPrepCubes import ff_utils
from OpenMMCubes import utils
import numpy as np
from oeommtools import utils as oeommutils


proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS',
                   'LEU', 'MET', 'PRO', 'THR', 'TYR',
                   'ARG', 'ASP', 'GLN', 'GLY', 'ILE',
                   'LYS', 'PHE', 'SER', 'TRP', 'VAL']

rnaResidues = ['A', 'G', 'C', 'U', 'I']
dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']


def applyffProtein(protein, opt):
    """
    This function applies the selected force field to the
    protein

    Parameters:
    -----------
    protein: OEMol molecule
        The protein to parametrize
    opt: python dictionary
        The options used to parametrize the protein

    Return:
    -------
    protein_structure: Parmed structure instance
        The parametrized protein parmed structure
    """

    topology, positions = oeommutils.oemol_to_openmmTop(protein)

    forcefield = app.ForceField(opt['protein_forcefield'])
    unmatched_residues = forcefield.getUnmatchedResidues(topology)

    if unmatched_residues:
        # Extended ff99SBildn force field
        oechem.OEThrow.Info("The following protein residues are not recognized "
                            "by the selected FF: {} - {}"
                            "\n...Extended FF is in use".format(opt['protein_forcefield'], unmatched_residues))

        ffext_fname = utils.get_data_filename('ComplexPrepCubes', 'ffext/amber99SBildn_ext.xml')
        forcefield = app.ForceField()
        forcefield.loadFile(ffext_fname)

        unmatched_residues = forcefield.getUnmatchedResidues(topology)

        if unmatched_residues:
            oechem.OEThrow.Fatal("Error. The following protein residues are not recognized "
                                 "by the extended force field {}".format(unmatched_residues))

    omm_system = forcefield.createSystem(topology, rigidWater=False)
    protein_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)

    return protein_structure


def applyffWater(water, opt):
    """
    This function applies the selected force field to the
    water

    Parameters:
    -----------
    water: OEMol molecule
        The water molecules to parametrize
    opt: python dictionary
        The options used to parametrize the water

    Return:
    -------
    water_structure: Parmed structure instance
        The parametrized water parmed structure
    """

    topology, positions = oeommutils.oemol_to_openmmTop(water)

    forcefield = app.ForceField(opt['solvent_forcefield'])
    unmatched_residues = forcefield.getUnmatchedResidues(topology)

    if unmatched_residues:
        oechem.OEThrow.Fatal("The following water molecules are not recognized "
                             "by the selected force field {}: {}".format(opt['solvent_forcefield'], unmatched_residues))

    omm_system = forcefield.createSystem(topology, rigidWater=False)
    water_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)

    return water_structure


def applyffExcipients(excipients, opt):
    """
    This function applies the selected force field to the
    excipients

    Parameters:
    -----------
    excipients: OEMol molecule
        The excipients molecules to parametrize
    opt: python dictionary
        The options used to parametrize the excipients

    Return:
    -------
    excipient_structure: Parmed structure instance
        The parametrized excipient parmed structure
    """

    # OpenMM topology and positions from OEMol
    topology, positions = oeommutils.oemol_to_openmmTop(excipients)

    # Try to apply the selected FF on the excipients
    forcefield = app.ForceField(opt['protein_forcefield'])

    # List of the unrecognized excipients
    unmatched_res_list = forcefield.getUnmatchedResidues(topology)

    # Unique unrecognized excipient names
    templates = set()
    for res in unmatched_res_list:
        templates.add(res.name)

    if templates:  # Some excipients are not recognized
        oechem.OEThrow.Info("The following excipients are not recognized "
                            "by the protein FF: {}"
                            "\nThey will be parametrized by using the FF: {}".format(templates, opt['other_forcefield']))

        # Create a bit vector mask used to split recognized from un-recognize excipients
        bv = oechem.OEBitVector(excipients.GetMaxAtomIdx())
        bv.NegateBits()

        # Dictionary containing the name and the parmed structures of the unrecognized excipients
        unrc_excipient_structures = {}

        # Dictionary used to skip already selected unrecognized excipients and count them
        unmatched_excp = {}

        # Ordered list of the unrecognized excipients
        unmatched_res_order = []

        for r_name in templates:
            unmatched_excp[r_name] = 0

        hv = oechem.OEHierView(excipients)

        for chain in hv.GetChains():
            for frag in chain.GetFragments():
                for hres in frag.GetResidues():
                    r_name = hres.GetOEResidue().GetName()
                    if r_name not in unmatched_excp:
                        continue
                    else:
                        unmatched_res_order.append(r_name)
                        if unmatched_excp[r_name]:  # Test if we have selected the unknown excipient
                            # Set Bit mask
                            atms = hres.GetAtoms()
                            for at in atms:
                                bv.SetBitOff(at.GetIdx())
                            unmatched_excp[r_name] += 1
                        else:
                            unmatched_excp[r_name] = 1
                            #  Create AtomBondSet to extract from the whole excipient system
                            #  the current selected FF unknown excipient
                            atms = hres.GetAtoms()
                            bond_set = set()
                            for at in atms:
                                bv.SetBitOff(at.GetIdx())
                                bonds = at.GetBonds()
                                for bond in bonds:
                                    bond_set.add(bond)
                            atom_bond_set = oechem.OEAtomBondSet(atms)
                            for bond in bond_set:
                                atom_bond_set.AddBond(bond)

                            # Create the unrecognized excipient OEMol
                            unrc_excp = oechem.OEMol()
                            if not oechem.OESubsetMol(unrc_excp, excipients, atom_bond_set):
                                oechem.OEThrow.Fatal("Is was not possible extract the residue: {}".format(r_name))

                            # Charge the unrecognized excipient
                            if not oequacpac.OEAssignCharges(unrc_excp,
                                                             oequacpac.OEAM1BCCCharges(symmetrize=True)):
                                oechem.OEThrow.Fatal("Is was not possible to "
                                                     "charge the extract residue: {}".format(r_name))

                            # If GAFF or GAFF2 is selected as FF check for tleap command
                            if opt['other_forcefield'] in ['GAFF', 'GAFF2']:
                                ff_utils.ParamLigStructure(oechem.OEMol(), opt['other_forcefield']).checkTleap

                            if opt['other_forcefield'] == 'SMIRNOFF':
                                unrc_excp = oeommutils.sanitizeOEMolecule(unrc_excp)

                            # Parametrize the unrecognized excipient by using the selected FF
                            pmd = ff_utils.ParamLigStructure(unrc_excp, opt['other_forcefield'],
                                                             prefix_name=opt['prefix_name']+'_'+r_name)
                            unrc_excp_struc = pmd.parameterize()
                            unrc_excp_struc.residues[0].name = r_name
                            unrc_excipient_structures[r_name] = unrc_excp_struc

        # Recognized FF excipients
        pred_rec = oechem.OEAtomIdxSelected(bv)
        rec_excp = oechem.OEMol()
        oechem.OESubsetMol(rec_excp, excipients, pred_rec)

        if rec_excp.NumAtoms() > 0:
            top_known, pos_known = oeommutils.oemol_to_openmmTop(rec_excp)
            ff_rec = app.ForceField(opt['protein_forcefield'])
            try:
                omm_system = ff_rec.createSystem(top_known, rigidWater=False)
                rec_struc = parmed.openmm.load_topology(top_known, omm_system, xyz=pos_known)
            except:
                oechem.OEThrow.Fatal("Error in the recognised excipient parametrization")

        # Unrecognized FF excipients
        bv.NegateBits()
        pred_unrc = oechem.OEAtomIdxSelected(bv)
        unrc_excp = oechem.OEMol()
        oechem.OESubsetMol(unrc_excp, excipients, pred_unrc)

        # Unrecognized FF excipients coordinates
        oe_coord_dic = unrc_excp.GetCoords()
        unrc_coords = np.ndarray(shape=(unrc_excp.NumAtoms(), 3))
        for at_idx in oe_coord_dic:
            unrc_coords[at_idx] = oe_coord_dic[at_idx]

        # It is important the order used to assemble the structures. In order to
        # avoid mismatch between the coordinates and the structures, it is convenient
        # to use the unrecognized residue order
        unmatched_res_order_count = []
        i = 0
        while i < len(unmatched_res_order):
            res_name = unmatched_res_order[i]
            for j in range(i+1, len(unmatched_res_order)):
                if unmatched_res_order[j] == res_name:
                    continue
                else:
                    break
            if i == (len(unmatched_res_order) - 1):
                num = 1
                unmatched_res_order_count.append((res_name, num))
                break
            else:
                num = j - i
                unmatched_res_order_count.append((res_name, num))
                i = j

        # Merge all the unrecognized Parmed structure
        unrc_struc = parmed.Structure()

        for pair in unmatched_res_order_count:
            res_name = pair[0]
            nums = pair[1]
            unrc_struc = unrc_struc + nums*unrc_excipient_structures[res_name]

        # Set the unrecognized coordinates
        unrc_struc.coordinates = unrc_coords

        # Set the parmed excipient structure merging
        # the unrecognized and recognized parmed
        # structures together
        if rec_excp.NumAtoms() > 0:
            excipients_structure = unrc_struc + rec_struc
        else:
            excipients_structure = unrc_struc

        return excipients_structure
    else:  # All the excipients are recognized by the selected FF
        omm_system = forcefield.createSystem(topology, rigidWater=False)
        excipients_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)

        return excipients_structure


# def applyffExcipients(excipients, opt):
#     """
#     This function applies the selected force field to the
#     excipients
#
#     Parameters:
#     -----------
#     excipients: OEMol molecule
#         The excipients molecules to parametrize
#     opt: python dictionary
#         The options used to parametrize the excipients
#
#     Return:
#     -------
#     excipient_structure: Parmed structure instance
#         The parametrized excipient parmed structure
#     """
#
#     # OpenMM topology and positions from OEMol
#     topology, positions = oeommutils.oemol_to_openmmTop(excipients)
#
#     # Try to apply the selected FF on the excipients
#     forcefield = app.ForceField(opt['protein_forcefield'])
#
#     # Unmatched residue templates
#     [templates, residues] = forcefield.generateTemplatesForUnmatchedResidues(topology)
#
#     if templates:  # Some excipients are not recognized
#         oechem.OEThrow.Info("The following excipients are not recognized "
#                             "by the protein FF: {}"
#                             "\nThey will be parametrized by using the FF: {}".format(residues, opt['other_forcefield']))
#
#         # Create a bit vector mask used to split recognized from un-recognize excipients
#         bv = oechem.OEBitVector(excipients.GetMaxAtomIdx())
#         bv.NegateBits()
#
#         # Dictionary containing the name and the parmed structures of the unrecognized excipients
#         unrc_excipient_structures = {}
#
#         # Dictionary used to skip already selected unrecognized excipients and count them
#         unmatched_excp = {}
#
#         # List to keep track of the unrecognized residue order
#         unmatched_res_order = []
#
#         for t in templates:
#             if t.name not in unmatched_excp:
#                 unmatched_excp[t.name] = 0
#
#         hv = oechem.OEHierView(excipients)
#
#         for chain in hv.GetChains():
#             for frag in chain.GetFragments():
#                 for hres in frag.GetResidues():
#                     r_name = hres.GetOEResidue().GetName()
#                     if r_name not in unmatched_excp:
#                         continue
#                     else:
#                         unmatched_res_order.append(r_name)
#                         if unmatched_excp[r_name]:  # Test if we have selected the unknown excipient
#                             # Set Bit mask
#                             atms = hres.GetAtoms()
#                             for at in atms:
#                                 bv.SetBitOff(at.GetIdx())
#                             unmatched_excp[r_name] += 1
#                         else:
#                             unmatched_excp[r_name] = 1
#                             #  Create AtomBondSet to extract from the whole excipient system
#                             #  the current selected FF unknown excipient
#                             atms = hres.GetAtoms()
#                             bond_set = set()
#                             for at in atms:
#                                 bv.SetBitOff(at.GetIdx())
#                                 bonds = at.GetBonds()
#                                 for bond in bonds:
#                                     bond_set.add(bond)
#                             atom_bond_set = oechem.OEAtomBondSet(atms)
#                             for bond in bond_set:
#                                 atom_bond_set.AddBond(bond)
#
#                             # Create the unrecognized excipient OEMol
#                             unrc_excp = oechem.OEMol()
#                             if not oechem.OESubsetMol(unrc_excp, excipients, atom_bond_set):
#                                 oechem.OEThrow.Fatal("Is was not possible extract the residue: {}".format(r_name))
#
#                             # Charge the unrecognized excipient
#                             if not oequacpac.OEAssignCharges(unrc_excp,
#                                                              oequacpac.OEAM1BCCCharges(symmetrize=True)):
#                                 oechem.OEThrow.Fatal("Is was not possible to "
#                                                      "charge the extract residue: {}".format(r_name))
#
#                             # If GAFF or GAFF2 is selected as FF check for tleap command
#                             if opt['other_forcefield'] in ['GAFF', 'GAFF2']:
#                                 ff_utils.ParamLigStructure(oechem.OEMol(), opt['other_forcefield']).checkTleap
#
#                             if opt['other_forcefield'] == 'SMIRNOFF':
#                                 unrc_excp = oeommutils.sanitizeOEMolecule(unrc_excp)
#                                 #oechem.OETriposAtomNames(unrc_excp)
#
#                             # Parametrize the unrecognized excipient by using the selected FF
#                             pmd = ff_utils.ParamLigStructure(unrc_excp, opt['other_forcefield'],
#                                                              prefix_name=opt['prefix_name']+'_'+r_name)
#                             unrc_excp_struc = pmd.parameterize()
#                             unrc_excp_struc.residues[0].name = r_name
#                             unrc_excipient_structures[r_name] = unrc_excp_struc
#
#         # Recognized FF excipients
#         pred_rec = oechem.OEAtomIdxSelected(bv)
#         rec_excp = oechem.OEMol()
#         oechem.OESubsetMol(rec_excp, excipients, pred_rec)
#         top_known, pos_known = oeommutils.oemol_to_openmmTop(rec_excp)
#         ff_rec = app.ForceField(opt['protein_forcefield'])
#
#         if rec_excp.NumAtoms() > 0:
#             try:
#                 omm_system = ff_rec.createSystem(top_known, rigidWater=False)
#                 rec_struc = parmed.openmm.load_topology(top_known, omm_system, xyz=pos_known)
#             except:
#                 oechem.OEThrow.Fatal("Error in the recognised excipient parametrization")
#
#         # Unrecognized FF excipients
#         bv.NegateBits()
#         pred_unrc = oechem.OEAtomIdxSelected(bv)
#         unrc_excp = oechem.OEMol()
#         oechem.OESubsetMol(unrc_excp, excipients, pred_unrc)
#
#         # ofs = oechem.oemolostream("gac.pdb")
#         # oechem.OEWriteConstMolecule(ofs, unrc_excp)
#
#         # Unrecognized FF excipients coordinates
#         oe_coord_dic = unrc_excp.GetCoords()
#         unrc_coords = np.ndarray(shape=(unrc_excp.NumAtoms(), 3))
#         for at_idx in oe_coord_dic:
#             unrc_coords[at_idx] = oe_coord_dic[at_idx]
#
#         # Merge all the unrecognized Parmed structure
#         unrc_struc = parmed.Structure()
#
#         # It is important the order used to assemble the structures. In order to
#         # avoid mismatch between the coordinates and the structures, it is convenient
#         # to use the unrecognized residue order
#
#         for res_name in unmatched_res_order:
#             unrc_struc = unrc_struc + unrc_excipient_structures[res_name]
#
#         # Set the unrecognized coordinates
#         unrc_struc.coordinates = unrc_coords
#
#         # Set the parmed excipient structure merging
#         # the unrecognized and recognized parmed
#         # structures together
#         if rec_excp.NumAtoms() > 0:
#             excipients_structure = unrc_struc + rec_struc
#         else:
#             excipients_structure = unrc_struc
#
#         return excipients_structure
#     else:  # All the excipients are recognized by the selected FF
#         omm_system = forcefield.createSystem(topology, rigidWater=False)
#         excipients_structure = parmed.openmm.load_topology(topology, omm_system, xyz=positions)
#
#         return excipients_structure


def applyffLigand(ligand, opt):
    """
    This function applies the selected force field to the
    ligand

    Parameters:
    -----------
    ligand: OEMol molecule
        The ligand molecule to parametrize
    opt: python dictionary
        The options used to parametrize the ligand

    Return:
    -------
    ligand_structure: Parmed structure instance
        The parametrized ligand parmed structure
    """

    # Check TLeap
    if opt['ligand_forcefield'] in ['GAFF', 'GAFF2']:
        ff_utils.ParamLigStructure(oechem.OEMol(), opt['ligand_forcefield']).checkTleap

    # Parametrize the Ligand
    pmd = ff_utils.ParamLigStructure(ligand, opt['ligand_forcefield'], prefix_name=opt['prefix_name'])
    ligand_structure = pmd.parameterize()
    ligand_structure.residues[0].name = "LIG"

    return ligand_structure


def order_check(mol, fname):
    """
    TO REMOVE
    This function is used to debug
    """
    import logging
    logger = logging.getLogger('Testing')
    hdlr = logging.FileHandler(fname)
    formatter = logging.Formatter('%(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)

    #hv = oechem.OEHierView(mol, oechem.OEAssumption_BondedResidue + oechem.OEAssumption_ResPerceived)
    hv = oechem.OEHierView(mol)

    for chain in hv.GetChains():
        logger.info('{}'.format(chain.GetChainID()))
        for frag in chain.GetFragments():
            for hres in frag.GetResidues():
                logger.info('\t{} {}'.format(hres.GetOEResidue().GetName(), hres.GetOEResidue().GetResidueNumber()))
                for oe_at in hres.GetAtoms():
                    logger.info('\t\t{} {}'.format(oe_at.GetName(), oe_at.GetIdx()))

    return
