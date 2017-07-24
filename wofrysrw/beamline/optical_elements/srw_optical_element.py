
from wofry.beamline.decorators import WOOpticalElementDecorator

from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters

from srwlib import SRWLOptC, srwl


class SRWOpticalElementDecorator:
    def toSRWLOpt(self):
        raise NotImplementedError("")

    @classmethod
    def fromSRWLOpt(cls, srwlopt=None):
        raise NotImplementedError("")

class SRWOpticalElement(SRWOpticalElementDecorator, WOOpticalElementDecorator):

    def applyOpticalElement(self, wavefront, parameters=None):

        if not parameters.has_additional_parameter("srw_oe_wavefront_propagation_parameters"):
            wavefront_propagation_parameters = WavefrontPropagationParameters()
        else:
            wavefront_propagation_parameters = parameters.get_additional_parameter("srw_oe_wavefront_propagation_parameters")

            if not isinstance(wavefront_propagation_parameters, WavefrontPropagationParameters):
                raise ValueError("SRW Wavefront Propagation Parameters not present")

        optBL = SRWLOptC([self.toSRWLOpt()],
                         [wavefront_propagation_parameters.to_SRW_array()])

        srwl.PropagElecField(wavefront, optBL)

        return wavefront

