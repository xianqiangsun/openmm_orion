# -*- coding: utf-8 -*-
import os

from tempfile import NamedTemporaryFile
import numpy as np
from floe.api.orion import in_orion, StreamingDataset

from simtk import unit
import yaml
from yank.analyze import get_analyzer
# Prevents repeated downloads of the same Dataset
download_cache = {}


def molecule_is_charged(mol):
    """
    Return True if the molecule has any nonzero charges; False otherwise.

    Parameters
    ----------
    mol : OEMol
       The molecule to examine.

    Returns
    -------
    is_charged : bool
        If any atom has nonzero charges, True is returned; else False.
    """
    is_charged = False
    for atom in mol.GetAtoms():
        if atom.GetPartialCharge() != 0.0:
            is_charged = True

    return is_charged
    

def get_data_filename(relative_path):
    """Get the full path to one of the reference files in testsystems.
    In the source distribution, these files are in ``examples/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.
    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).
    """

    from pkg_resources import resource_filename
    fn = resource_filename('examples', os.path.join('data', relative_path))
    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)
    return fn


def download_dataset_to_file(dataset_id):
    """
    Used to retrieve a dataset either from Orion or from the local machine
    """
    if in_orion():
        if dataset_id in download_cache:
            return download_cache[dataset_id]
        if os.path.isfile(dataset_id):
            download_cache[dataset_id] = dataset_id
            return dataset_id
        tmp = NamedTemporaryFile(suffix=".oeb.gz", delete=False)
        stream = StreamingDataset(dataset_id, input_format=".oeb.gz")
        stream.download_to_file(tmp.name)
        download_cache[dataset_id] = tmp.name
        return tmp.name
    else:
        return dataset_id


def unfuck_oechem_mol2_file(filename, substructure_name='MOL'):
    """Undo oechem fuckery with mol2 substructure names.

    Parameters
    ----------
    filename : str
        mol2 filename for file to be modified
    substructure_name : str, optional, default='MOL'
        substructure name to replace `<0>`

    """
    # Replace <0> substructure names with valid text.
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    newlines = [line.replace('<0>', substructure_name) for line in lines]
    outfile = open(filename, 'w')
    outfile.writelines(newlines)
    outfile.close()


def analyze_directory(source_directory):
    """
    This Function has been copied and adapted from the Yank ver 0.17.0 source code
    (yank.analyse.analyze_directory)

    Analyze contents of store files to compute free energy differences.

    This function is needed to preserve the old auto-analysis style of YANK. What it exactly does can be refined
    when more analyzers and simulations are made available. For now this function exposes the API.

    Parameters
    ----------
    source_directory : string
        The location of the simulation storage files.

    """
    analysis_script_path = os.path.join(source_directory, 'analysis.yaml')
    if not os.path.isfile(analysis_script_path):
        err_msg = 'Cannot find analysis.yaml script in {}'.format(source_directory)
        raise RuntimeError(err_msg)
    with open(analysis_script_path, 'r') as f:
        analysis = yaml.load(f)

    data = dict()
    for phase_name, sign in analysis:
        phase_path = os.path.join(source_directory, phase_name + '.nc')
        phase = get_analyzer(phase_path)
        data[phase_name] = phase.analyze_phase()
        kT = phase.kT

    # Compute free energy and enthalpy
    DeltaF = 0.0
    dDeltaF = 0.0
    DeltaH = 0.0
    dDeltaH = 0.0
    for phase_name, sign in analysis:
        DeltaF -= sign * (data[phase_name]['DeltaF'] + data[phase_name]['DeltaF_standard_state_correction'])
        dDeltaF += data[phase_name]['dDeltaF'] ** 2
        DeltaH -= sign * (data[phase_name]['DeltaH'] + data[phase_name]['DeltaF_standard_state_correction'])
        dDeltaH += data[phase_name]['dDeltaH'] ** 2
    dDeltaF = np.sqrt(dDeltaF)
    dDeltaH = np.sqrt(dDeltaH)

    DeltaF = DeltaF * kT / unit.kilocalories_per_mole
    dDeltaF = dDeltaF * kT / unit.kilocalories_per_mole
    DeltaH = DeltaH * kT / unit.kilocalories_per_mole
    dDeltaH = dDeltaH * kT / unit.kilocalories_per_mole

    return DeltaF, dDeltaF, DeltaH, dDeltaH