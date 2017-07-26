import numpy

from syned.beamline.optical_elements.mirrors.mirror import Mirror
from syned.beamline.shape import Ellipsoid, Rectangle, Plane

from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement
from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters

from srwlib import SRWLOptC, SRWLOptMir
from srwlib import srwl, srwl_opt_setup_surf_height_1d, srwl_opt_setup_surf_height_2d, srwl_uti_read_data_cols

class Orientation:
    UP = 0
    DOWN = 1
    LEFT = 2
    RIGHT = 3

class ApertureShape:
    RECTANGULAR = 'r'
    ELLIPTIC = 'e'

class SimulationMethod:
    THIN = 1
    THICK = 2

class TreatInputOutput:
    WAVEFRONT_INPUT_PLANE_BEFORE_OUTPUT_PLANE_AFTER  = 0
    WAVEFRONT_INPUT_CENTER_OUTPUT_CENTER  = 1
    WAVEFRONT_INPUT_CENTER_OUTPUT_CENTER_DRIFT_BACK  = 2

class ScaleType:
    LINEAR = 'lin'
    LOGARITHMIC = 'log'

class SRWMirror(Mirror, SRWOpticalElement):
    def __init__(self,
                 name                               = "Undefined",
                 tangential_size                    = 1.2,
                 sagittal_size                      = 0.01,
                 grazing_angle                      = 0.003,
                 orientation_of_reflection_plane    = Orientation.UP,
                 height_profile_data_file           = "mirror.dat",
                 height_profile_data_file_dimension = 1,
                 height_amplification_coefficient   = 1.0):


        Mirror.__init__(self,
                        name=name,
                        boundary_shape=Rectangle(x_left=-0.5*sagittal_size,
                                                 x_right=0.5*sagittal_size,
                                                 y_bottom=-0.5*tangential_size,
                                                 y_top=0.5*tangential_size),
                        surface_shape=Plane())

        self.tangential_size                                  = tangential_size
        self.sagittal_size                                    = sagittal_size
        self.grazing_angle                                    = grazing_angle
        self.orientation_of_reflection_plane                  = orientation_of_reflection_plane

        self.height_profile_data_file = height_profile_data_file
        self.height_profile_data_file_dimension = height_profile_data_file_dimension
        self.height_amplification_coefficient = height_amplification_coefficient

    def applyOpticalElement(self, wavefront, parameters=None):
        wavefront = super().applyOpticalElement(wavefront, parameters)

        if not self.height_profile_data_file is None:

            if self.height_profile_data_file_dimension == 1:
                height_profile_data = srwl_uti_read_data_cols(self.height_profile_data_file,
                                                              _str_sep='\t',
                                                              _i_col_start=0,
                                                              _i_col_end=1)

                optTrEr = srwl_opt_setup_surf_height_1d(_height_prof_data=height_profile_data,
                                                        _ang=self.grazing_angle,
                                                        _dim='y',
                                                        _amp_coef=self.height_amplification_coefficient)
            elif self.height_profile_data_file_dimension == 2:
                height_profile_data = srwl_uti_read_data_cols(self.height_profile_data_file,
                                                              _str_sep='\t',
                                                              _i_col_start=0,
                                                              _i_col_end=2)

                optTrEr = srwl_opt_setup_surf_height_2d(_height_prof_data=height_profile_data,
                                                        _ang=self.grazing_angle,
                                                        _dim='y',
                                                        _amp_coef=self.height_amplification_coefficient)

            optBL = SRWLOptC([optTrEr],
                             [WavefrontPropagationParameters().to_SRW_array()])

            srwl.PropagElecField(wavefront, optBL)

        return wavefront

    def toSRWLOpt(self):
        if self.orientation_of_reflection_plane == Orientation.LEFT:
            nvx = numpy.cos(self.grazing_angle)
            nvy = 0
            nvz = -numpy.sin(self.grazing_angle)
            tvx = -numpy.sin(self.grazing_angle)
            tvy = 0
        elif self.orientation_of_reflection_plane == Orientation.RIGHT:
            nvx = numpy.cos(self.grazing_angle)
            nvy = 0
            nvz = -numpy.sin(self.grazing_angle)
            tvx = numpy.sin(self.grazing_angle)
            tvy = 0
        elif self.orientation_of_reflection_plane == Orientation.UP:
            nvx = 0
            nvy = numpy.cos(self.grazing_angle)
            nvz = -numpy.sin(self.grazing_angle)
            tvx = 0
            tvy = numpy.sin(self.grazing_angle)
        elif self.orientation_of_reflection_plane == Orientation.DOWN:
            nvx = 0
            nvy = numpy.cos(self.grazing_angle)
            nvz = -numpy.sin(self.grazing_angle)
            tvx = 0
            tvy = -numpy.sin(self.grazing_angle)

        return self.get_SRWLOptMir(nvx, nvy, nvz, tvx, tvy)

    def get_SRWLOptMir(self, nvx, nvy, nvz, tvx, tvy):
        return SRWLOptMir(_size_tang=self.tangential_size,
                          _size_sag=self.sagittal_size,
                          _ap_shape=ApertureShape.RECTANGULAR,
                          _sim_meth=SimulationMethod.THICK,
                          _treat_in_out=TreatInputOutput.WAVEFRONT_INPUT_CENTER_OUTPUT_CENTER,
                          _nvx=nvx,
                          _nvy=nvy,
                          _nvz=nvz,
                          _tvx=tvx,
                          _tvy=tvy)

    def fromSRWLOpt(self, srwlopt=SRWLOptMir()):
        if not isinstance(srwlopt, SRWLOptMir):
            raise ValueError("SRW object is not a SRWLOptMir object")

        if srwlopt.tvx != 0.0:
            orientation_of_reflection_plane = Orientation.LEFT if srwlopt.tvx < 0 else Orientation.RIGHT
            grazing_angle = abs(numpy.arctan(srwlopt.nvz/srwlopt.nvx))
        else:
            orientation_of_reflection_plane = Orientation.DOWN if srwlopt.tvy < 0 else Orientation.UP
            grazing_angle = abs(numpy.arctan(srwlopt.nvz/srwlopt.nvy))

        self.__init__(tangential_size                 = srwlopt.dt,
                      sagittal_size                   = srwlopt.ds,
                      grazing_angle                   = grazing_angle,
                      orientation_of_reflection_plane = orientation_of_reflection_plane,
                      height_profile_data_file        = None)

