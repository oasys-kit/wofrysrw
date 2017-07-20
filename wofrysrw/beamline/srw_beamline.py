
from syned.beamline.beamline import Beamline

from wofrysrw.storage_ring.srw_light_source import SRWLightSource
from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters

class SRWBeamline(Beamline):

    def __init__(self,
                 light_source=SRWLightSource(),
                 beamline_elements_list=[],
                 wavefront_propagation_parameters_list=[]):
        super().__init__(light_source=light_source, beamline_elements_list=beamline_elements_list)

        self._wavefront_propagation_parameters_list = wavefront_propagation_parameters_list

    # overwrites the SynedObject method for dealing with list
    def to_dictionary(self):
        dict_to_save = super().to_dictionary()
        dict_to_save["srw_wavefront_propagation_parameters"] = [ el.to_dictionary() for el in self._wavefront_propagation_parameters_list ]

        return dict_to_save

    def append_wavefront_propagation_parameters(self, wavefront_propagation_parameters=WavefrontPropagationParameters()):
        if not isinstance(wavefront_propagation_parameters,WavefrontPropagationParameters):
            raise Exception("Input class must be of type: "+WavefrontPropagationParameters.__name__)
        else:
            self._wavefront_propagation_parameters_list.append(wavefront_propagation_parameters)

    def get_wavefront_propagation_parameters(self):
        return len(self._wavefront_propagation_parameters_list)

    def get_wavefront_propagation_parameters_at(self, index):
        if index >= len(self._wavefront_propagation_parameters_list):
            raise IndexError("Index " + str(index) + " out of bounds")

        return self._wavefront_propagation_parameters_list[index]
