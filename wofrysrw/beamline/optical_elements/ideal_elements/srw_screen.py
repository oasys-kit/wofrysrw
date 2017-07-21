from syned.beamline.optical_elements.ideal_elements.screen import Screen

from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement

class SRWScreen(Screen, SRWOpticalElement):
    def __init__(self, name="Undefined"):
        Screen.__init__(self, name=name)

    def applyOpticalElement(self, wavefront, parameters=None):
        return wavefront
