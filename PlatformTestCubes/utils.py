# -*- coding: utf-8 -*-
import os

from tempfile import NamedTemporaryFile

from openeye.oechem import(
    oemolostream, OEWriteConstMolecule
)
from openeye.oedocking import OEWriteReceptorFile

from floe.api.orion import in_orion, StreamingDataset

# Prevents repeated downloads of the same Dataset
download_cache = {}


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


def dump_query(prefix, name, qmol, receptor):
    """
    Writes the Molecule or receptor out to file on the machine
    """
    tag = "{0}_{1}.query".format(prefix, name)
    query_file = "{0}.oeb.gz".format(tag)
    with oemolostream(query_file) as ofs:
        OEWriteConstMolecule(ofs, qmol)
    if receptor.IsValid():
        receptor_file = "{0}.receptor.oeb.gz".format(tag)
        OEWriteReceptorFile(receptor, receptor_file)
    return tag, query_file
