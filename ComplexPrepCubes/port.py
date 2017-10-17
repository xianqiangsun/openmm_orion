from floe.api import (parameter, MoleculeOutputPort, SourceCube, OEMolOStreamCube)
from floe.api.orion import StreamingDataset, config_from_env, MultipartDatasetUploader
from openeye import oechem
# from floe.api.ports import BinaryMoleculeInputPort, MoleculeSerializerMixin
# from floe.api.parameter import DataSetOutputParameter


class ProteinReader(SourceCube):
    title = "Protein Reader Cube"
    version = "0.0.0"
    classification = [["Protein Reader Cube", "OEChem", "Reader Cube"]]
    tags = ['OEChem']
    description = """
    A Protein Reader Cube 
    Input:
    -------
    oechem.OEMCMol or - Streamed-in of the protein system
    The input file can be an .oeb, .oeb.gz, .pdb or a .mol2 file

    Output:
    -------
    oechem.OEMCMol - Emits the protein system
    """

    success = MoleculeOutputPort("success")

    data_in = parameter.DataSetInputParameter(
        "data_in",
        help_text="Protein to read in",
        required=True,
        description="The Protein to read in")

    limit = parameter.IntegerParameter(
        "limit",
        required=False)

    download_format = parameter.StringParameter(
        "download_format",
        choices=[".oeb.gz", ".oeb", ".pdb", ".mol2", ".smi"],
        required=False,
        default=".oeb.gz")

    protein_prefix = parameter.StringParameter(
        'protein_prefix',
        default='PRT',
        help_text='The protein prefix name used to identify the protein')

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
                    mol.SetTitle(self.opt['protein_prefix'])
                    yield mol
                    count += 1
                    if max_idx is not None and count == max_idx:
                        break
        else:
            stream = StreamingDataset(self.args.data_in,
                                      input_format=self.args.download_format)
            for mol in stream:
                mol.SetTitle(self.opt['protein_prefix'])
                yield mol
                count += 1
                if max_idx is not None and count == max_idx:
                    break


# class SimOutputCube(OEMolOStreamCube):
#     """
#     A sink cube that writes molecules to a file
#     """
#     classification = [["Output"]]
#     title = "Output Writer"
#
#     intake = BinaryMoleculeInputPort('intake')
#     data_out = DataSetOutputParameter('data_out',
#                                       required=True,
#                                       title='Name of Dataset to create',
#                                       description='The dataset to output')
#     backend = DataSetOutputParameter(
#         'backend',
#         default="auto",
#         choices=["db", "s3", "auto"],
#         description="The Orion storage backend to use")
#
#     def begin(self):
#         self.in_orion = config_from_env() is not None
#         self.decoder = MoleculeSerializerMixin()
#         self.need_decode = not self.args.data_out.endswith(".oeb.gz")
#         if self.in_orion:
#             self.ofs = MultipartDatasetUploader(self.args.data_out,
#                                                 tags=[self.name],
#                                                 backend=self.args.backend)
#         elif self.need_decode:
#             self.ofs = oechem.oemolostream(str(self.args.data_out))
#         else:
#             self.ofs = open(str(self.args.data_out), 'wb')
#
#     def write(self, mol, port):
#         if self.in_orion or not self.need_decode:
#             self.ofs.write(mol)
#         else:
#             oechem.OEWriteMolecule(self.ofs, self.decoder.decode(mol))
#
#     def end(self):
#         if self.in_orion:
#             self.ofs.complete()
#         else:
#             self.ofs.close()