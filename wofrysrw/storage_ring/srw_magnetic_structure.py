from syned.storage_ring.magnetic_structure import MagneticStructure
from srwlib import array, SRWLMagFldC

from wofrysrw.srw_object import SRWObject

class SRWMagneticStructureDecorator():

    def get_SRWMagneticStructure(self):
        raise NotImplementedError("this method should be implented in subclasses")

    def get_SRWLMagFldC(self):
        return SRWLMagFldC([self.get_SRWMagneticStructure()], array('d', [0]), array('d', [0]), array('d', [0]))


class SRWMagneticStructure(SRWMagneticStructureDecorator, SRWObject):
    def __init__(self,
                 horizontal_central_position = 0.0,
                 vertical_central_position = 0.0,
                 longitudinal_central_position = 0.0):
        super().__init__()

        self.horizontal_central_position = horizontal_central_position
        self.vertical_central_position = vertical_central_position
        self.longitudinal_central_position = longitudinal_central_position

    def to_python_code(self, data=None):
        text_code  = self.to_python_code_aux()
        text_code += "magnetic_field_container = SRWLMagFldC([magnetic_structure], " + \
                     "array('d', [" + str(self.horizontal_central_position) + "]), " + \
                     "array('d', [" + str(self.vertical_central_position) + "])), " + \
                     "array('d', [" + str(self.longitudinal_central_position) + "])))"  + "\n"

        return text_code

    def to_python_code_aux(self):
        raise NotImplementedError()
