
from wofry.beamline.decorators import OpticalElementDecorator

from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters, WavefrontPropagationOptionalParameters

from srwlib import SRWLOptC, srwl


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

