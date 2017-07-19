
import numpy

from syned.beamline.optical_elements.absorbers.slit import Slit
from syned.beamline.shape import BoundaryShape, Rectangle, Ellipse

from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement

from srwlib import SRWLOptA


class SRWSlit(Slit, SRWOpticalElement):
    def __init__(self, name="Undefined", boundary_shape=BoundaryShape()):
        Slit.__init__(self, name=name, boundary_shape=boundary_shape)

    def applySRWOpticalElement(self, wavefront, parameters=None):
        boundaries = self._boundary_shape.get_boundaries()

        if isinstance(self._boundary_shape, Rectangle):
            wavefront.mesh.xStart = boundaries[0] #Initial Horizontal Position [m]
            wavefront.mesh.xFin = boundaries[1] #Final Horizontal Position [m]
            wavefront.mesh.yStart =  boundaries[2] #Initial Vertical Position [m]
            wavefront.mesh.yFin = boundaries[3] #Final Vertical Position [m]
        elif isinstance(self._boundary_shape, Ellipse):
            wavefront.mesh.xStart = boundaries[0] #Initial Horizontal Position [m]
            wavefront.mesh.xFin = boundaries[1] #Final Horizontal Position [m]
            wavefront.mesh.yStart =  boundaries[2] #Initial Vertical Position [m]
            wavefront.mesh.yFin = boundaries[3] #Final Vertical Position [m]
        else:
            raise ValueError("Wrong Boundary Shape type")

    def toSRWLOpt(self):
        boundaries = self._boundary_shape.get_boundaries()

        Dx = numpy.abs(boundaries[1]-boundaries[0])
        Dy = numpy.abs(boundaries[3]-boundaries[2])
        x = 0.5*(boundaries[1]-boundaries[0])
        y = 0.5*(boundaries[3]-boundaries[2])

        if isinstance(self._boundary_shape, Rectangle):
            shape = 'r'
        elif isinstance(self._boundary_shape, Ellipse):
            if Dx != Dy:
                raise ValueError("SRW doesn't support elliptic apertures")

            shape = 'c'

        return SRWLOptA(_shape=shape,
                        _ap_or_ob='a',
                        _Dx=Dx,
                        _Dy=Dy,
                        _x=x,
                        _y=y)


    @classmethod
    def fromSRWLOpt(cls, srwlopt=SRWLOptA()):
        if not isinstance(srwlopt, SRWLOptA):
            raise ValueError("SRW object is not a SRWLOptA object")
        
        if not srwlopt.ap_or_ob == 'a':
            raise ValueError("SRW object is an obstruction")

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

        return Slit(boundary_shape=boundary_shape)