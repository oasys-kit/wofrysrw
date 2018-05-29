import numpy

from syned.beamline.optical_elements.mirrors.mirror import Mirror
from syned.beamline.shape import Rectangle, Ellipse, Circle

from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement, Orientation
from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters

from srwlib import SRWLOptC, SRWLOptMir
from srwlib import srwl, srwl_opt_setup_surf_height_1d, srwl_opt_setup_surf_height_2d, srwl_uti_read_data_cols


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
                 vertical_position_of_mirror_center = 0.0,
                 horizontal_position_of_mirror_center = 0.0,
                 grazing_angle                      = 0.003,
                 orientation_of_reflection_plane    = Orientation.UP,
                 invert_tangent_component           = False,
                 height_profile_data_file           = "mirror.dat",
                 height_profile_data_file_dimension = 1,
                 height_amplification_coefficient   = 1.0):


        Mirror.__init__(self,
                        name=name,
                        boundary_shape=Rectangle(x_left=horizontal_position_of_mirror_center - 0.5*sagittal_size,
                                                 x_right=horizontal_position_of_mirror_center + 0.5*sagittal_size,
                                                 y_bottom=vertical_position_of_mirror_center - 0.5*tangential_size,
                                                 y_top=vertical_position_of_mirror_center + 0.5*tangential_size),
                        surface_shape=self.get_shape())

        self.tangential_size                                  = tangential_size
        self.sagittal_size                                    = sagittal_size
        self.grazing_angle                                    = grazing_angle
        self.orientation_of_reflection_plane                  = orientation_of_reflection_plane
        self.invert_tangent_component                         = invert_tangent_component

        self.height_profile_data_file = height_profile_data_file
        self.height_profile_data_file_dimension = height_profile_data_file_dimension
        self.height_amplification_coefficient = height_amplification_coefficient

    def get_shape(self):
        raise NotImplementedError()

    def applyOpticalElement(self, wavefront, parameters=None):
        wavefront = super().applyOpticalElement(wavefront, parameters)

        if not self.height_profile_data_file is None:

            if self.orientation_of_reflection_plane == Orientation.LEFT or self.orientation_of_reflection_plane == Orientation.RIGHT:
                dim = 'x'
            elif self.orientation_of_reflection_plane == Orientation.UP or self.orientation_of_reflection_plane == Orientation.DOWN:
                dim = 'y'

            if self.height_profile_data_file_dimension == 1:
                height_profile_data = srwl_uti_read_data_cols(self.height_profile_data_file,
                                                              _str_sep='\t',
                                                              _i_col_start=0,
                                                              _i_col_end=1)

                optTrEr = srwl_opt_setup_surf_height_1d(_height_prof_data=height_profile_data,
                                                        _ang=self.grazing_angle,
                                                        _dim=dim,
                                                        _amp_coef=self.height_amplification_coefficient)
            elif self.height_profile_data_file_dimension == 2:
                height_profile_data = srwl_uti_read_data_cols(self.height_profile_data_file,
                                                              _str_sep='\t')

                optTrEr = srwl_opt_setup_surf_height_2d(_height_prof_data=height_profile_data,
                                                        _ang=self.grazing_angle,
                                                        _dim=dim,
                                                        _amp_coef=self.height_amplification_coefficient)

            optBL = SRWLOptC([optTrEr],
                             [WavefrontPropagationParameters().to_SRW_array()])

            srwl.PropagElecField(wavefront, optBL)

        return wavefront

    def toSRWLOpt(self):
        nvx, nvy, nvz, tvx, tvy = self.get_orientation_vectors()
        x, y = self.getXY()

        if isinstance(self.get_boundary_shape(), Rectangle):
            ap_shape = ApertureShape.RECTANGULAR
        elif isinstance(self.get_boundary_shape(), Ellipse) or isinstance(self.get_boundary_shape(), Circle):
            ap_shape = ApertureShape.ELLIPTIC

        return self.get_SRWLOptMir(nvx, nvy, nvz, tvx, tvy, x, y, ap_shape)

    def get_SRWLOptMir(self, nvx, nvy, nvz, tvx, tvy, x, y, ap_shape):
        mirror = SRWLOptMir()

        mirror.set_dim_sim_meth(_size_tang=self.tangential_size,
                                _size_sag=self.sagittal_size,
                                _ap_shape=ap_shape,
                                _sim_meth=SimulationMethod.THICK,
                                _treat_in_out=TreatInputOutput.WAVEFRONT_INPUT_CENTER_OUTPUT_CENTER)
        mirror.set_orient(_nvx=nvx,
                          _nvy=nvy,
                          _nvz=nvz,
                          _tvx=tvx,
                          _tvy=tvy,
                          _x = x,
                          _y = y)

        return mirror

    def fromSRWLOpt(self, srwlopt=SRWLOptMir()):
        if not isinstance(srwlopt, SRWLOptMir):
            raise ValueError("SRW object is not a SRWLOptMir object")

        if srwlopt.tvx != 0.0:
            orientation_of_reflection_plane = Orientation.LEFT if srwlopt.nvx < 0 else Orientation.RIGHT
            grazing_angle = abs(numpy.arctan(srwlopt.nvz/srwlopt.nvx))
            if orientation_of_reflection_plane == Orientation.LEFT: invert_tangent_component = numpy.sign(srwlopt.nvz) == numpy.sign(srwlopt.tvx)
            else: invert_tangent_component = numpy.sign(srwlopt.nvz) != numpy.sign(srwlopt.tvx)
        elif srwlopt.tvy != 0.0:
            orientation_of_reflection_plane = Orientation.DOWN if srwlopt.nvy < 0 else Orientation.UP
            grazing_angle = abs(numpy.arctan(srwlopt.nvz/srwlopt.nvy))
            if orientation_of_reflection_plane == Orientation.UP: invert_tangent_component = numpy.sign(srwlopt.nvy) == numpy.sign(srwlopt.tvy)
            else: invert_tangent_component = numpy.sign(srwlopt.nvy) != numpy.sign(srwlopt.tvy)
        else:
            raise ValueError("Tangential orientation angles (tvx/tvy) are both 0.0!")

        self.__init__(tangential_size                 = srwlopt.dt,
                      sagittal_size                   = srwlopt.ds,
                      grazing_angle                   = grazing_angle,
                      orientation_of_reflection_plane = orientation_of_reflection_plane,
                      invert_tangent_component        = invert_tangent_component,
                      height_profile_data_file        = None)
