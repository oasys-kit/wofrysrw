

from syned.beamline.beamline import Beamline

from wofrysrw.storage_ring.srw_light_source import SRWLightSource


#Meaning of Wavefront Propagation Parameters:
#[0]: Auto-Resize (1) or not (0) Before propagation
#[1]: Auto-Resize (1) or not (0) After propagation
#[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
#[3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
#[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
#[5]: Horizontal Range modification factor at Resizing (1. means no modification)
#[6]: Horizontal Resolution modification factor at Resizing
#[7]: Vertical Range modification factor at Resizing
#[8]: Vertical Resolution modification factor at Resizing
#[9]: Type of wavefront Shift before Resizing (not yet implemented)
#[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
#[11]: New Vertical wavefront Center position after Shift (not yet implemented)

class SRWWavefrontPropagationParameters(object):
    def __init__(self,
                 auto_resize_before_propagation                         = 0,
                 auto_resize_after_propagation                          = 0,
                 relative_prevision_for_propagation_with_autoresizing   = 1.0,
                 allow_semianalytical_treatment_of_quadratic_phase_term = 0,
                 do_any_resizing_on_fourier_side_using_fft              = 0,
                 horizontal_range_modification_factor_at_resizing       = 1.0,
                 horizontal_resolution_modification_factor_at_resizing  = 1.0,
                 vertical_range_modification_factor_at_resizing         = 1.0,
                 vertical_resolution_modification_factor_at_resizing    = 1.0,
                 type_of_wavefront_shift_before_resizing                = 0,
                 new_horizontal_wavefront_center_position_after_shift   = 0,
                 new_vertical_wavefront_center_position_after_shift     = 0):
            self._auto_resize_before_propagation                         = auto_resize_before_propagation,
            self._auto_resize_after_propagation                          = auto_resize_after_propagation,
            self._relative_prevision_for_propagation_with_autoresizing   = relative_prevision_for_propagation_with_autoresizing,
            self._allow_semianalytical_treatment_of_quadratic_phase_term = allow_semianalytical_treatment_of_quadratic_phase_term,
            self._do_any_resizing_on_fourier_side_using_fft              = do_any_resizing_on_fourier_side_using_fft,
            self._horizontal_range_modification_factor_at_resizing       = horizontal_range_modification_factor_at_resizing,
            self._horizontal_resolution_modification_factor_at_resizing  = horizontal_resolution_modification_factor_at_resizing,
            self._vertical_range_modification_factor_at_resizing         = vertical_range_modification_factor_at_resizing,
            self._vertical_resolution_modification_factor_at_resizing    = vertical_resolution_modification_factor_at_resizing
            self._type_of_wavefront_shift_before_resizing                = type_of_wavefront_shift_before_resizing,
            self._new_horizontal_wavefront_center_position_after_shift   = new_horizontal_wavefront_center_position_after_shift,
            self._new_vertical_wavefront_center_position_after_shift     = new_vertical_wavefront_center_position_after_shift

    def to_SRW_array(self):
        return [self._auto_resize_before_propagation,
                self._auto_resize_after_propagation,
                self._relative_prevision_for_propagation_with_autoresizing,
                self._allow_semianalytical_treatment_of_quadratic_phase_term,
                self._do_any_resizing_on_fourier_side_using_fft,
                self._horizontal_range_modification_factor_at_resizing,
                self._horizontal_resolution_modification_factor_at_resizing,
                self._vertical_range_modification_factor_at_resizing,
                self._vertical_resolution_modification_factor_at_resizing,
                self._type_of_wavefront_shift_before_resizing,
                self._new_horizontal_wavefront_center_position_after_shift,
                self._new_vertical_wavefront_center_position_after_shift]

class SRWBeamline(Beamline):

    def __init__(self,
                 light_source=SRWLightSource(),
                 beamline_elements_list=[],
                 srw_wavefront_propagation_parameters_list=[]):
        super().__init__(light_source=light_source, beamline_elements_list=beamline_elements_list)

        self._srw_wavefront_propagation_parameters_list = srw_wavefront_propagation_parameters_list


    # overwrites the SynedObject method for dealing with list
    def to_dictionary(self):
        dict_to_save = super().to_dictionary()
        dict_to_save["srw_wavefront_propagation_parameters"] = [ el.to_dictionary() for el in self._srw_wavefront_propagation_parameters_list ]

        return dict_to_save

    def append_srw_wavefront_propagation_parameters(self, srw_wavefront_propagation_parameters=SRWWavefrontPropagationParameters()):
        if not isinstance(srw_wavefront_propagation_parameters,SRWWavefrontPropagationParameters):
            raise Exception("Input class must be of type: "+SRWWavefrontPropagationParameters.__name__)
        else:
            self._srw_wavefront_propagation_parameters_list.append(srw_wavefront_propagation_parameters)

    def get_srw_wavefront_propagation_parameters(self):
        return len(self._srw_wavefront_propagation_parameters_list)

    def get_srw_wavefront_propagation_parameters_at(self, index):
        if index >= len(self._srw_wavefront_propagation_parameters_list):
            raise IndexError("Index " + str(index) + " out of bounds")

        return self._srw_wavefront_propagation_parameters_list[index]
