from srwlib import SRWLMagFldM

from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet
from wofrysrw.storage_ring.srw_magnetic_structure import SRWMagneticStructureDecorator

class SRWBendingMagnet(BendingMagnet, SRWMagneticStructureDecorator):

    def __init__(self,
                 radius = 0.0,
                 magnetic_field = 0.0,
                 length = 0.0):
        BendingMagnet.__init__(self, radius, magnetic_field, length)

    def get_SRWMagneticStructure(self):
        return SRWLMagFldM(self._magnetic_field, 1, 'n', self._length)


