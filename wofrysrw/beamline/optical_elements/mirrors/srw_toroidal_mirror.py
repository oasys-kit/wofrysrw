from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement
from wofrysrw.beamline.optical_elements.mirrors.srw_mirror import SRWMirror, Orientation, ApertureShape, SimulationMethod, TreatInputOutput

from srwlib import SRWLOptMirTor

class SRWToroidalMirror(SRWMirror, SRWOpticalElement):
    def __init__(self,
                 name                                       = "Undefined",
                 tangential_size                            = 1.2,
                 sagittal_size                              = 0.01,
                 grazing_angle                              = 0.003,
                 orientation_of_reflection_plane            = Orientation.UP,
                 tangential_radius =1,
                 sagittal_radius=1,
                 height_profile_data_file                   = "mirror.dat",
                 height_profile_data_file_dimension         = 1,
                 height_amplification_coefficient           = 1.0):

        super().__init__(name=name,
                         tangential_size=tangential_size,
                         sagittal_size=sagittal_size,
                         grazing_angle=grazing_angle,
                         orientation_of_reflection_plane=orientation_of_reflection_plane,
                         height_profile_data_file=height_profile_data_file,
                         height_profile_data_file_dimension=height_profile_data_file_dimension,
                         height_amplification_coefficient=height_amplification_coefficient)

        self.tangential_radius  = tangential_radius
        self.sagittal_radius = sagittal_radius


    def get_SRWLOptMir(self, nvx, nvy, nvz, tvx, tvy):
        return SRWLOptMirTor(_size_tang=self.tangential_size,
                            _size_sag=self.sagittal_size,
                            _rt=self.tangential_radius,
                            _rs=self.sagittal_radius,
                            _ang_graz=self.grazing_angle,
                            _ap_shape=ApertureShape.RECTANGULAR,
                            _sim_meth=SimulationMethod.THICK,
                            _treat_in_out=TreatInputOutput.WAVEFRONT_INPUT_CENTER_OUTPUT_CENTER,
                            _nvx=nvx,
                            _nvy=nvy,
                            _nvz=nvz,
                            _tvx=tvx,
                            _tvy=tvy)

    def fromSRWLOpt(self, srwlopt=SRWLOptMirTor()):
        if not isinstance(srwlopt, SRWLOptMirTor):
            raise ValueError("SRW object is not a SRWLOptMirEl object")

        super().fromSRWLOpt(srwlopt)

        self.tangential_radius = srwlopt.radTan
        self.sagittal_radius = srwlopt.radSag
