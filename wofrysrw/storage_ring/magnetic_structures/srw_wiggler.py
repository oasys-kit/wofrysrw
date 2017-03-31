from srwlib import array, SRWLMagFld3D, SRWLMagFldH, SRWLMagFldU

from syned.storage_ring.magnetic_structures.wiggler import Wiggler
from wofrysrw.storage_ring.srw_magnetic_structure import SRWMagneticStructureDecorator

import numpy

class SRWWiggler(Wiggler, SRWMagneticStructureDecorator):

    def __init__(self,
                 K_vertical=0.0,
                 K_horizontal=0.0,
                 period_length = 0.0,
                 number_of_periods = 1):
        Wiggler.__init__(self,
                         K_vertical=K_vertical,
                         K_horizontal=K_horizontal,
                         period_length=period_length,
                         number_of_periods=number_of_periods)

    def get_SRWMagneticStructure(self):
        if self._number_of_periods <=3:
            n_points = 100

            arBx = array('d', [0]*n_points) if self._K_horizontal > 0.0 else None
            arBy = array('d', [0]*n_points) if self._K_vertical > 0.0 else None

            longitudinal_mesh = numpy.linspace(-self._period_length/2., self._period_length/2., n_points)
            phases = 2.0*numpy.pi*(longitudinal_mesh/self._period_length)

            for i in range (0, n_points):
                if self._K_horizontal > 0.0: arBx[i] = -self.magnetic_field_horizontal() * numpy.cos(phases[i])
                if self._K_vertical > 0.0:   arBy[i] = self.magnetic_field_vertical() * numpy.cos(phases[i])

            return SRWLMagFld3D(_arBx=arBx,
                                _arBy=arBy,
                                _nx=1,
                                _ny=1,
                                _nz=n_points,
                                _rz=self._period_length,
                                _nRep=int(self._number_of_periods),
                                _interp=1)
        else:
            magnetic_fields = []

            if self._K_vertical > 0.0:
                magnetic_fields.append(SRWLMagFldH(1, 'v', self.magnetic_field_vertical(), 0, 1, 1))

            if self._K_horizontal > 0.0:
                magnetic_fields.append(SRWLMagFldH(1, 'h', self.magnetic_field_horizontal(), 0, -1, 1))

            return SRWLMagFldU(magnetic_fields,
                               self._period_length,
                               int(self._number_of_periods))