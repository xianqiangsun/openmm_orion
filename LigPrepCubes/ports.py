from __future__ import unicode_literals
from floe.api import OutputPort, InputPort, Port
from floe.constants import BYTES
from floe.api import (parameter, MoleculeOutputPort, SourceCube)
from floe.api.orion import StreamingDataset, config_from_env

from openeye import oechem

try:
    import cPickle as pickle
except ImportError:
    import pickle


class MoleculeSerializerMixin(object):

    def __init__(self, *args, **kwargs):
        super(MoleculeSerializerMixin, self).__init__(*args, **kwargs)
        self._ifs = oechem.oemolistream()
        self._ifs.SetFormat(oechem.OEFormat_OEB)
        self._ifs.Setgz(True)
        errs = oechem.oeosstream()
        self._ofs = oechem.oemolostream(errs, False)
        self._ofs.openstring()
        self._ofs.Setgz(True)
        self._ofs.SetFormat(oechem.OEFormat_OEB)

    def encode(self, *mols):
        """
        By default, serializes molecules as gzipped oeb files in a string
        """
        for mol in mols:
            code = oechem.OEWriteMolecule(self._ofs, mol)
            if code != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Unable to encode mol: {}".format(code))
        res = self._ofs.GetString()
        self._ofs.close()
        self._ofs.openstring()
        return res

    def decode(self, mol_data):
        """
        By default, deserializes data into molecules for use in the cube
        """
        mol = oechem.OEMol()
        if type(mol_data) == oechem.OEMol:
            return mol_data
        if not self._ifs.openstring(mol_data):
            raise RuntimeError("Failed to open string")
        if not oechem.OEReadMolecule(self._ifs, mol):
            print("Unable to decode molecule")
        self._ifs.close()
        return mol


class MoleculePortSerializer(MoleculeSerializerMixin, Port):

    OEB_GZ = '.oeb.gz'
    PORT_TYPES = (BYTES, OEB_GZ)
    FORMAT = OEB_GZ


class CustomMoleculeInputPort(InputPort, MoleculePortSerializer):
    pass


class CustomMoleculeOutputPort(OutputPort, MoleculePortSerializer):
    pass


class LigandReader(SourceCube):
    title = "LigandReader Cube"
    version = "0.0.0"
    classification = [["Ligand Reader Cube", "OEChem", "Reader Cube"]]
    tags = ['OEChem']
    description = """
    Ligand Reader Cube 
    Input:
    -------
    oechem.OEMCMol or - Streamed-in of Ligands
    The input file can be an .oeb, .oeb.gz, .pdb or a .mol2 file

    Output:
    -------
    oechem.OEMCMol - Emits the Ligands
    """

    success = MoleculeOutputPort("success")

    data_in = parameter.DataSetInputParameter(
        "data_in",
        help_text="Ligand to read in",
        required=True,
        description="The Ligand to read in")

    limit = parameter.IntegerParameter(
        "limit",
        required=False)

    download_format = parameter.StringParameter(
        "download_format",
        choices=[".oeb.gz", ".oeb", ".pdb", ".mol2", ".smi"],
        required=False,
        default=".oeb.gz")

    prefix = parameter.StringParameter(
        'prefix',
        default='',
        help_text='An SD tag used as prefix string')

    suffix = parameter.StringParameter(
        'suffix',
        default='',
        help_text='An SD tag used as suffix string')

    type = parameter.StringParameter(
        'type',
        default='LIG',
        required=True,
        help_text='The ligand reside name')

    IDTag = parameter.BooleanParameter(
        'IDTag',
        default=True,
        required=False,
        help_text='If True/Checked ligands are enumerated by sequentially integers.'
                  'A SD tag containing part of the ligand name and an integer is used '
                  'to create a unique IDTag which is attached to the ligand')

    def begin(self):
        self.opt = vars(self.args)

    def __iter__(self):
        max_idx = self.args.limit
        if max_idx is not None:
            max_idx = int(max_idx)
        count = 0
        self.config = config_from_env()
        in_orion = self.config is not None
        if not in_orion:
            with oechem.oemolistream(str(self.args.data_in)) as ifs:
                for mol in ifs.GetOEMols():
                    mol.SetData(oechem.OEGetTag('prefix'), self.opt['prefix'])
                    mol.SetData(oechem.OEGetTag('suffix'), self.opt['suffix'])

                    for at in mol.GetAtoms():
                        residue = oechem.OEAtomGetResidue(at)
                        residue.SetName(self.opt['type'])
                        oechem.OEAtomSetResidue(at, residue)

                    if self.opt['IDTag']:
                        mol.SetData(oechem.OEGetTag('IDTag'), 'l_' + mol.GetTitle()[0:12] + str(count))
                    yield mol
                    count += 1
                    if max_idx is not None and count == max_idx:
                        break
        else:
            stream = StreamingDataset(self.args.data_in,
                                      input_format=self.args.download_format)
            for mol in stream:
                mol.SetData(oechem.OEGetTag('prefix'), self.opt['prefix'])
                mol.SetData(oechem.OEGetTag('suffix'), self.opt['suffix'])

                for at in mol.GetAtoms():
                    residue = oechem.OEAtomGetResidue(at)
                    residue.SetName(self.opt['type'])
                    oechem.OEAtomSetResidue(at, residue)

                if self.opt['IDTag']:
                    mol.SetData(oechem.OEGetTag('IDTag'), 'l_' + mol.GetTitle()[0:12] + str(count))
                yield mol
                count += 1
                if max_idx is not None and count == max_idx:
                    break