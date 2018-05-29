from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement
from wofrysrw.beamline.optical_elements.mirrors.srw_mirror import SRWMirror, Orientation, ApertureShape, SimulationMethod, TreatInputOutput
from syned.beamline.shape import Plane

from srwlib import SRWLOptMirPl

class SRWPlaneMirror(SRWMirror, SRWOpticalElement):
    def __init__(self,
                 name                               = "Undefined",
                 tangential_size                    = 1.2,
                 sagittal_size                      = 0.01,
                 grazing_angle                      = 0.003,
                 orientation_of_reflection_plane    = Orientation.UP,
                 invert_tangent_component           = False,
                 height_profile_data_file           = "mirror.dat",
                 height_profile_data_file_dimension = 1,
                 height_amplification_coefficient   = 1.0):

        super().__init__(name=name,
                         tangential_size=tangential_size,
                         sagittal_size=sagittal_size,
                         grazing_angle=grazing_angle,
                         orientation_of_reflection_plane=orientation_of_reflection_plane,
                         invert_tangent_component=invert_tangent_component,
                         height_profile_data_file=height_profile_data_file,
                         height_profile_data_file_dimension=height_profile_data_file_dimension,
                         height_amplification_coefficient=height_amplification_coefficient)

    def get_shape(self):
        return Plane()

    def get_SRWLOptMir(self, nvx, nvy, nvz, tvx, tvy, x, y, ap_shape):
        return SRWLOptMirPl(_size_tang=self.tangential_size,
                            _size_sag=self.sagittal_size,
                            _ap_shape=ap_shape,
                            _sim_meth=SimulationMethod.THICK,
                            _treat_in_out=TreatInputOutput.WAVEFRONT_INPUT_CENTER_OUTPUT_CENTER,
                            _nvx=nvx,
                            _nvy=nvy,
                            _nvz=nvz,
                            _tvx=tvx,
                            _tvy=tvy,
                            _x=x,
                            _y=y)

    def fromSRWLOpt(self, srwlopt=SRWLOptMirPl()):
        if not isinstance(srwlopt, SRWLOptMirPl):
            raise ValueError("SRW object is not a SRWLOptMirEl object")

        super().fromSRWLOpt(srwlopt)
