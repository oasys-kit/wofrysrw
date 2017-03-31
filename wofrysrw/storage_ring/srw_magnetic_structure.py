from syned.storage_ring.magnetic_structure import MagneticStructure
from srwlib import array, SRWLMagFldC

class SRWMagneticStructureDecorator():

    def get_SRWMagneticStructure(self):
        raise NotImplementedError("this method should be implented in subclasses")

    def get_SRWLMagFldC(self):
        return SRWLMagFldC([self.get_SRWMagneticStructure()], array('d', [0]), array('d', [0]), array('d', [0]))

class SRWMagneticStructure(MagneticStructure, SRWMagneticStructureDecorator):
    def __init__(self):
        super().__init__()

