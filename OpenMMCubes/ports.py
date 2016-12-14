from __future__ import unicode_literals
"""
"""
import struct
from floe.api import OutputPort, InputPort
from floe.constants import BYTES
from simtk import openmm, unit

class OpenMMSystemSerializationMixin(object):

    def encode(self, openmm_system):
        return openmm.XmlSerializer.serialize(openmm_system).encode()

    def decode(self, serialized_system):
        return openmm.XmlSerializer.deserialize(serialized_system.decode())

class OpenMMSystemInput(OpenMMSystemSerializationMixin, InputPort):
    PORT_TYPES = (BYTES, "OpenMMSystem")

class OpenMMSystemOutput(OpenMMSystemSerializationMixin, OutputPort):
    PORT_TYPES = (BYTES, "OpenMMSystem")
