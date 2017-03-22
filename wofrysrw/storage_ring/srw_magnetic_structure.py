

from syned.storage_ring.magnetic_structure import MagneticStructure

class SRWMagneticStructureDecorator():

    def get_SRWLMagFldC(self):
        pass


class SRWMagneticStructure(MagneticStructure, SRWMagneticStructureDecorator):
    def __init__(self):
        super().__init__()
