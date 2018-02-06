import numpy

from wofry.beamline.decorators import OpticalElementDecorator

from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters, WavefrontPropagationOptionalParameters

from srwlib import SRWLOptC, srwl

class Orientation:
    UP = 0
    DOWN = 1
    LEFT = 2
    RIGHT = 3

class SRWOpticalElementDecorator:
    def toSRWLOpt(self):
        raise NotImplementedError("")

    @classmethod
    def fromSRWLOpt(cls, srwlopt=None):
        raise NotImplementedError("")

class SRWOpticalElement(SRWOpticalElementDecorator, OpticalElementDecorator):

    def applyOpticalElement(self, wavefront=None, parameters=None):

        if not parameters.has_additional_parameter("srw_oe_wavefront_propagation_parameters"):
            wavefront_propagation_parameters = WavefrontPropagationParameters()
        else:
            wavefront_propagation_parameters = parameters.get_additional_parameter("srw_oe_wavefront_propagation_parameters")

            if not isinstance(wavefront_propagation_parameters, WavefrontPropagationParameters):
                raise ValueError("SRW Wavefront Propagation Parameters are inconsistent")

        srw_parameters_array = wavefront_propagation_parameters.to_SRW_array()

        if parameters.has_additional_parameter("srw_oe_wavefront_propagation_optional_parameters"):
            wavefront_propagation_optional_parameters = parameters.get_additional_parameter("srw_oe_wavefront_propagation_optional_parameters")

            if not isinstance(wavefront_propagation_parameters, WavefrontPropagationOptionalParameters):
                raise ValueError("SRW Wavefront Propagation Optional Parameters are inconsistent")

            wavefront_propagation_optional_parameters.append_to_srw_array(srw_parameters_array)

        optBL = SRWLOptC([self.toSRWLOpt()],
                         [srw_parameters_array])

        srwl.PropagElecField(wavefront, optBL)

        return wavefront

    def get_orientation_vectors(self):
        if self.orientation_of_reflection_plane == Orientation.LEFT:
            nvx = numpy.cos(self.grazing_angle)
            nvy = 0
            nvz = -numpy.sin(self.grazing_angle)
            tvx = numpy.sin(self.grazing_angle)
            tvy = 0
        elif self.orientation_of_reflection_plane == Orientation.RIGHT:
            nvx = numpy.cos(self.grazing_angle)
            nvy = 0
            nvz = -numpy.sin(self.grazing_angle)
            tvx = -numpy.sin(self.grazing_angle)
            tvy = 0
        elif self.orientation_of_reflection_plane == Orientation.UP:
            nvx = 0
            nvy = numpy.cos(self.grazing_angle)
            nvz = -numpy.sin(self.grazing_angle)
            tvx = 0
            tvy = numpy.sin(self.grazing_angle)
        elif self.orientation_of_reflection_plane == Orientation.DOWN:
            nvx = 0
            nvy = -numpy.cos(self.grazing_angle)
            nvz = numpy.sin(self.grazing_angle)
            tvx = 0
            tvy = numpy.sin(self.grazing_angle)

        return nvx, nvy, nvz, tvx, tvy