# TODO: REMOVE THIS!!!!
try:
    from srwlib import *
    SRWLIB_AVAILABLE = True
except:
    try:
        from wpg.srwlib import *
        SRWLIB_AVAILABLE = True
    except:
        SRWLIB_AVAILABLE = False
        print("SRW is not available")

import scipy.constants as codata
angstroms_to_eV = codata.h*codata.c/codata.e*1e10

from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
from wofry.propagator.propagator import Propagator2D
from wofry.propagator.test.srw_wavefront import WOSRWWavefront

from wofrysrw.beamline.srw_beamline import SRWWavefrontPropagationParameters

class FresnelSRW(Propagator2D):

    HANDLER_NAME = "FRESNEL_SRW"

    def get_handler_name(self):
        return self.HANDLER_NAME

    """
    2D Fresnel propagator using convolution via Fourier transform
    :param wavefront:
    :param propagation_distance:
    :param srw_autosetting:set to 1 for automatic SRW redimensionate wavefront
    :return:
    """

    def do_specific_progation(self, wavefront, propagation_distance, parameters):
        if not SRWLIB_AVAILABLE: raise ImportError("SRW is not available")

        if not parameters.has_additional_parameter("srw_drift_wavefront_propagation_parameters"):
            srw_wavefront_propagation_parameters = SRWWavefrontPropagationParameters()
        else:
            srw_wavefront_propagation_parameters = parameters.get_additional_parameter("srw_drift_wavefront_propagation_parameters")

            if not isinstance(srw_wavefront_propagation_parameters, SRWWavefrontPropagationParameters):
                raise ValueError("SRW Wavefront Propagation Parameters not present")

        is_generic_wavefront = isinstance(wavefront, GenericWavefront2D)

        if is_generic_wavefront:
            wavefront = WOSRWWavefront.fromGenericWavefront(wavefront)
        else:
            if not isinstance(wavefront, WOSRWWavefront): raise ValueError("wavefront cannot be managed by this propagator")

        #
        # propagation (simple wavefront drift
        #

        optBL = SRWLOptC([SRWLOptD(propagation_distance)], # drift space
                         [srw_wavefront_propagation_parameters.to_SRW_array()]) #"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)

        srwl.PropagElecField(wavefront, optBL)

        if is_generic_wavefront:
            return wavefront.toGenericWavefront()
        else:
            return wavefront