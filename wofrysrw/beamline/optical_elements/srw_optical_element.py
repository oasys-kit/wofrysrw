import numpy

from syned.beamline.shape import Ellipse, Rectangle, Circle
from wofry.beamline.decorators import OpticalElementDecorator

from wofry.propagator.propagator import PropagationParameters
from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters, WavefrontPropagationOptionalParameters
from wofrysrw.srw_object import SRWObject

from srwlib import SRWLOptC, srwl

class Orientation:
    UP = 0
    DOWN = 1
    LEFT = 2
    RIGHT = 3

class SRWOpticalElementDecorator(SRWObject):
    def toSRWLOpt(self):
        raise NotImplementedError("")

    def fromSRWLOpt(self, srwlopt=None):
        raise NotImplementedError("")

class SRWOpticalElement(SRWOpticalElementDecorator, OpticalElementDecorator):

    def applyOpticalElement(self, wavefront=None, parameters=None):
        optBL = SRWLOptC([self.toSRWLOpt()],
                         [self.get_srw_wavefront_propagation_parameter(parameters)])

        srwl.PropagElecField(wavefront, optBL)

        return wavefront

    def add_to_srw_native_array(self, oe_array = [], pp_array=[], parameters=None, wavefront=None):
        oe_array.append(self.toSRWLOpt())
        pp_array.append(self.get_srw_wavefront_propagation_parameter(parameters))

    def get_srw_wavefront_propagation_parameter(self, parameters):
        if isinstance(parameters, PropagationParameters):
            if not parameters.has_additional_parameter("srw_oe_wavefront_propagation_parameters"):
                wavefront_propagation_parameters = WavefrontPropagationParameters()
            else:
                wavefront_propagation_parameters = parameters.get_additional_parameter("srw_oe_wavefront_propagation_parameters")

                if not isinstance(wavefront_propagation_parameters, WavefrontPropagationParameters):
                    raise ValueError("SRW Wavefront Propagation Parameters are inconsistent")

            srw_parameters_array = wavefront_propagation_parameters.to_SRW_array()

            if parameters.has_additional_parameter("srw_oe_wavefront_propagation_optional_parameters"):
                wavefront_propagation_optional_parameters = parameters.get_additional_parameter("srw_oe_wavefront_propagation_optional_parameters")

                if not isinstance(wavefront_propagation_optional_parameters, WavefrontPropagationOptionalParameters):
                    raise ValueError("SRW Wavefront Propagation Optional Parameters are inconsistent")

                wavefront_propagation_optional_parameters.append_to_srw_array(srw_parameters_array)
        elif isinstance(parameters, list):
            wavefront_propagation_parameters = parameters[0]
            wavefront_propagation_optional_parameters = parameters[1]

            srw_parameters_array = wavefront_propagation_parameters.to_SRW_array()
            if not wavefront_propagation_optional_parameters is None: wavefront_propagation_optional_parameters.append_to_srw_array(srw_parameters_array)

        return srw_parameters_array

    def getXY(self):
        boundary_shape = self.get_boundary_shape()

        if not boundary_shape is None:
            if isinstance(boundary_shape, Rectangle) or isinstance(boundary_shape, Ellipse):
                x_left, x_right, y_bottom, y_top = boundary_shape.get_boundaries()

                return x_left + 0.5*(x_right-x_left), \
                       y_bottom + 0.5*(y_top-y_bottom)

            elif isinstance(boundary_shape, Circle):
                radius, x_center, y_center = boundary_shape.get_boundaries()

                return x_center, y_center
        else:
            return 0.0, 0.0

    def get_orientation_vectors(self):
        sign = (-1 if self.invert_tangent_component else 1)

        if self.orientation_of_reflection_plane == Orientation.LEFT:
            nvx = -numpy.cos(self.grazing_angle)
            nvy = 0
            nvz = -numpy.sin(self.grazing_angle)
            tvx = sign*numpy.sin(self.grazing_angle)
            tvy = 0
        elif self.orientation_of_reflection_plane == Orientation.RIGHT:
            nvx = numpy.cos(self.grazing_angle)
            nvy = 0
            nvz = -numpy.sin(self.grazing_angle)
            tvx = -sign*numpy.sin(self.grazing_angle)
            tvy = 0
        elif self.orientation_of_reflection_plane == Orientation.UP:
            nvx = 0
            nvy = numpy.cos(self.grazing_angle)
            nvz = -numpy.sin(self.grazing_angle)
            tvx = 0
            tvy = -sign*numpy.sin(self.grazing_angle)
        elif self.orientation_of_reflection_plane == Orientation.DOWN:
            nvx = 0
            nvy = -numpy.cos(self.grazing_angle)
            nvz = -numpy.sin(self.grazing_angle)
            tvx = 0
            tvy = -sign*numpy.sin(self.grazing_angle)

        return nvx, nvy, nvz, tvx, tvy
