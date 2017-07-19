
import numpy

from syned.beamline.optical_elements.absorbers.slit import Slit
from syned.beamline.shape import BoundaryShape, Rectangle, Ellipse

from wofry.beamline.decorators import WOOpticalElementDecorator

from wofrysrw.beamline.srw_beamline import SRWWavefrontPropagationParameters

import srwlib as srwl
from srwlib import SRWLOptA, SRWLOptC


class SRWOpticalElement(WOOpticalElementDecorator):

    def applyOpticalElement(self, wavefront, parameters=None):

        self.applySRWOpticalElement(wavefront, parameters)

        if not parameters.has_additional_parameter("srw_oe_wavefront_propagation_parameters"):
            srw_wavefront_propagation_parameters = SRWWavefrontPropagationParameters()
        else:
            srw_wavefront_propagation_parameters = parameters.get_additional_parameter("srw_oe_wavefront_propagation_parameters")

            if not isinstance(srw_wavefront_propagation_parameters, SRWWavefrontPropagationParameters):
                raise ValueError("SRW Wavefront Propagation Parameters not present")

        optBL = SRWLOptC([self.toSRWLOpt()],
                         [srw_wavefront_propagation_parameters.to_SRW_array()]) #"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)

        srwl.PropagElecField(wavefront, optBL)

        return wavefront

    def applySRWOpticalElement(self, wavefront, parameters=None):
        raise NotImplementedError("")

    def toSRWLOpt(self):
        raise NotImplementedError("")

    @classmethod
    def fromSRWLOpt(cls, srwlopt=None):
        raise NotImplementedError("")