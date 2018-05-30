
from syned.beamline.beamline import Beamline

from wofrysrw.storage_ring.srw_light_source import SRWLightSource
from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters, WavefrontPropagationOptionalParameters

from srwlib import *

class Where:
    DRIFT_BEFORE = "before"
    DRIFT_AFTER = "after"
    OE = "oe"

    @classmethod
    def tuple(cls):
        return cls.DRIFT_BEFORE, cls.OE, cls.DRIFT_AFTER

class SRWBeamline(Beamline):

    def __init__(self,
                 light_source=SRWLightSource(),
                 beamline_elements_list=[]):
        super().__init__(light_source=light_source, beamline_elements_list=beamline_elements_list)

        self._wavefront_propagation_parameters_list = {}
        self._wavefront_propagation_parameters_list[Where.DRIFT_BEFORE] = []
        self._wavefront_propagation_parameters_list[Where.OE]           = []
        self._wavefront_propagation_parameters_list[Where.DRIFT_AFTER]  = []

    # overwrites the SynedObject method for dealing with list
    def to_dictionary(self): # TODO: to be adjusted....
        dict_to_save = super().to_dictionary()
        dict_to_save["srw_wavefront_propagation_parameters"] = [ el[0].to_dictionary() for el in self._wavefront_propagation_parameters_list[Where.OE]]

        return dict_to_save

    def append_wavefront_propagation_parameters(self, wavefront_propagation_parameters=WavefrontPropagationParameters(), wavefront_propagation_optional_parameters=WavefrontPropagationOptionalParameters(), where=Where.OE):
        self._wavefront_propagation_parameters_list[where].append([wavefront_propagation_parameters, wavefront_propagation_optional_parameters])

    def get_wavefront_propagation_parameters(self, where=Where.OE):
        return self._wavefront_propagation_parameters_list[where]

    def get_wavefront_propagation_parameters_at(self, index, where=Where.OE):
        if index >= len(self._wavefront_propagation_parameters_list[where]):
            raise IndexError("Index " + str(index) + " out of bounds")

        return self._wavefront_propagation_parameters_list[where][index]

    def duplicate(self):
        beamline_elements_list = []
        for beamline_element in self._beamline_elements_list:
            beamline_elements_list.append(beamline_element)

        new_beamline = SRWBeamline(light_source=self._light_source, beamline_elements_list = beamline_elements_list)

        for where in Where.tuple():
            for element in self.get_wavefront_propagation_parameters(where):
                new_beamline.append_wavefront_propagation_parameters(element[0], element[1], where)

        return new_beamline

    def to_python_code(self):
        return ""
