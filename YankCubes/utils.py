# -*- coding: utf-8 -*-
import os

from tempfile import NamedTemporaryFile

from openeye.oechem import(
    oemolostream, OEWriteConstMolecule
)
from openeye import oechem
from openeye.oedocking import OEWriteReceptorFile
import numpy as np
from floe.api.orion import in_orion, StreamingDataset
from simtk.openmm.app import Topology
from simtk.openmm.app.element import Element
from simtk import unit, openmm
# Prevents repeated downloads of the same Dataset
download_cache = {}

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
