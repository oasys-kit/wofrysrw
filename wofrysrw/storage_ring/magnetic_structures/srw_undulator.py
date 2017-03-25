from srwlib import array, SRWLMagFldH, SRWLMagFldU, SRWLMagFldC

from syned.storage_ring.magnetic_structures.undulator import Undulator
from wofrysrw.storage_ring.srw_magnetic_structure import SRWMagneticStructureDecorator

class SRWUndulator(Undulator, SRWMagneticStructureDecorator):

    def __init__(self,
                 K_vertical = 0.0,
                 K_horizontal = 0.0,
                 period_length = 0.0,
                 number_of_periods = 1):
        Undulator.__init__(self, K_vertical, K_horizontal, period_length, number_of_periods)

    def get_SRWLMagFldC(self):

        magnetic_fields = []

        if self._K_vertical > 0.0:
            magnetic_fields.append(SRWLMagFldH(1, 'v', self.B_vertical(), 0, 1, 1))

        if self._K_horizontal > 0.0:
            magnetic_fields.append(SRWLMagFldH(1, 'h', self.B_horizontal(), 0, -1, 1))

        srw_undulator = SRWLMagFldU(magnetic_fields,
                                    self._period_length,
                                    self._number_of_periods)

        return SRWLMagFldC([srw_undulator], array('d', [0]), array('d', [0]), array('d', [0]))
