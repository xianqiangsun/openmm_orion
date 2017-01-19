from __future__ import unicode_literals
from floe.api import OutputPort, InputPort
from floe.constants import BYTES
from simtk import openmm, unit
from lzma import compress, decompress



class OpenMMSystemSerializationMixin(object):

    def encode(self, openmm_system):
        serialized_system = openmm.XmlSerializer.serialize(openmm_system)
        encoded = serialized_system.encode()
        #return openmm.XmlSerializer.serialize(openmm_system).encode()
        return compress( encoded )

    def decode(self, compressed_system):
        decompresed = decompress(compressed_system)
        decoded = decompresed.decode()
        system = openmm.XmlSerializer.deserialize( decoded )
        return system
        #return decompress( openmm.XmlSerializer.deserialize(serialized_system).decode() )


class OpenMMSystemInput(OpenMMSystemSerializationMixin, InputPort):
    PORT_TYPES = (BYTES, "OpenMMSystem")

class OpenMMSystemOutput(OpenMMSystemSerializationMixin, OutputPort):
    PORT_TYPES = (BYTES, "OpenMMSystem")
