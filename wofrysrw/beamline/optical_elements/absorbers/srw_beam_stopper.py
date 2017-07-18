
import numpy

from syned.beamline.optical_elements.absorbers.beam_stopper import BeamStopper
from syned.beamline.shape import BoundaryShape, Rectangle, Ellipse

from wofry.beamline.decorators import WOOpticalElementDecorator

from srwlib import SRWLOptA

class SRWBeamStopper(BeamStopper, WOOpticalElementDecorator):
    def __init__(self, name="Undefined", boundary_shape=BoundaryShape()):
        BeamStopper.__init__(self, name=name, boundary_shape=boundary_shape)

    def applyOpticalElement(self, wavefront):
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

        return wavefront

    def toSRWLOptA(self):
        boundaries = self._boundary_shape.get_boundaries()

        Dx = numpy.abs(boundaries[1]-boundaries[0])
        Dy = numpy.abs(boundaries[3]-boundaries[2])
        x = 0.5*(boundaries[1]-boundaries[0])
        y = 0.5*(boundaries[3]-boundaries[2])

        if isinstance(self._boundary_shape, Rectangle):
            shape = 'r'
        elif isinstance(self._boundary_shape, Ellipse):
            if Dx != Dy:
                raise ValueError("SRW doesn't support elliptic obstructions")

            shape = 'c'

        return SRWLOptA(_shape=shape,
                        _ap_or_ob='o',
                        _Dx=Dx,
                        _Dy=Dy,
                        _x=x,
                        _y=y)


    @classmethod
    def fromSRWLOptA(cls, srwlopta=SRWLOptA()):
        if not isinstance(srwlopta, SRWLOptA):
            raise ValueError("SRW object is not a SRWLOptA object")

        if not srwlopta.ap_or_ob == 'o':
            raise ValueError("SRW object is an aperture")

        if srwlopta.shape == 'r':
            boundary_shape = Rectangle(x_left=srwlopta.x - 0.5*srwlopta.Dx,
                                       x_right=srwlopta.x + 0.5*srwlopta.Dx,
                                       y_bottom=srwlopta.y - 0.5*srwlopta.Dy,
                                       y_top=srwlopta.y + 0.5*srwlopta.Dy)
        elif srwlopta.shape == 'c':
            boundary_shape = Ellipse(min_ax_left=srwlopta.x - 0.5*srwlopta.Dx,
                                     min_ax_right=srwlopta.x + 0.5*srwlopta.Dx,
                                     maj_ax_bottom=srwlopta.y - 0.5*srwlopta.Dy,
                                     maj_ax_top=srwlopta.y + 0.5*srwlopta.Dy)

        return BeamStopper(boundary_shape=boundary_shape)

