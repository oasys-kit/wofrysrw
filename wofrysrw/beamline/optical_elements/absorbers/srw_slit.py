from syned.beamline.optical_elements.absorbers.slit import Slit
from syned.beamline.shape import BoundaryShape, Rectangle, Ellipse

from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement

from srwlib import SRWLOptA

class SRWSlit(Slit, SRWOpticalElement):
    def __init__(self, name="Undefined", boundary_shape=BoundaryShape()):
        Slit.__init__(self, name=name, boundary_shape=boundary_shape)

    def toSRWLOpt(self):
        boundaries = self._boundary_shape.get_boundaries()

        Dx = abs(boundaries[1]-boundaries[0])
        Dy = abs(boundaries[3]-boundaries[2])
        x = 0.5*(boundaries[1]+boundaries[0])
        y = 0.5*(boundaries[3]+boundaries[2])

        if isinstance(self._boundary_shape, Rectangle):
            shape = 'r'
        elif isinstance(self._boundary_shape, Ellipse):
            if Dx != Dy:
                raise ValueError("SRW doesn't support elliptic apertures")

            shape = 'c'

        return SRWLOptA(_shape=shape,
                        _ap_or_ob=self.get_srw_ap_or_ob(),
                        _Dx=Dx,
                        _Dy=Dy,
                        _x=x,
                        _y=y)


    def fromSRWLOpt(self, srwlopt=SRWLOptA()):
        if not isinstance(srwlopt, SRWLOptA):
            raise ValueError("SRW object is not a SRWLOptA object")
        
        if not srwlopt.ap_or_ob == self.get_srw_ap_or_ob():
            raise ValueError("SRW object not compatible")

        if srwlopt.shape == 'r':
            boundary_shape = Rectangle(x_left=srwlopt.x - 0.5 * srwlopt.Dx,
                                       x_right=srwlopt.x + 0.5 * srwlopt.Dx,
                                       y_bottom=srwlopt.y - 0.5 * srwlopt.Dy,
                                       y_top=srwlopt.y + 0.5 * srwlopt.Dy)
        elif srwlopt.shape == 'c':
            boundary_shape = Ellipse(min_ax_left=srwlopt.x - 0.5 * srwlopt.Dx,
                                     min_ax_right=srwlopt.x + 0.5 * srwlopt.Dx,
                                     maj_ax_bottom=srwlopt.y - 0.5 * srwlopt.Dy,
                                     maj_ax_top=srwlopt.y + 0.5 * srwlopt.Dy)

        self.__init__(boundary_shape=boundary_shape)

    def get_srw_ap_or_ob(self):
        raise NotImplementedError()