from syned.storage_ring.magnetic_structure import MagneticStructure
from srwlib import array, SRWLMagFldC

from wofrysrw.srw_object import SRWObject

class SRWMagneticStructureDecorator():

    def get_SRWMagneticStructure(self):
        raise NotImplementedError("this method should be implented in subclasses")

    def get_SRWLMagFldC(self):
        return SRWLMagFldC([self.get_SRWMagneticStructure()], array('d', [0]), array('d', [0]), array('d', [0]))


class SRWMagneticStructure(SRWMagneticStructureDecorator, SRWObject):
    def __init__(self):
        super().__init__()

    def to_python_code(self, data=None):
        text_code  = self.to_python_code_aux()
        text_code += "magnetic_field_container = SRWLMagFldC([magnetic_structure], array('d', [0]), array('d', [0]), array('d', [0]))"  + "\n"

        return text_code

    def to_python_code_aux(self):
        raise NotImplementedError()
