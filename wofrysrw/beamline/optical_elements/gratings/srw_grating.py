import numpy

from syned.beamline.optical_elements.gratings.grating import Grating
from syned.beamline.shape import Ellipse, Rectangle, Circle

from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement
from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters

from srwlib import SRWLOptC, SRWLOptMir, SRWLOptG
from srwlib import srwl, srwl_opt_setup_surf_height_1d, srwl_opt_setup_surf_height_2d, srwl_uti_read_data_cols

from wofrysrw.beamline.optical_elements.mirrors.srw_mirror import Orientation, TreatInputOutput, ApertureShape, SimulationMethod

class SRWGrating(Grating, SRWOpticalElement):
    def __init__(self,
                 name                               = "Undefined",
                 tangential_size                    = 1.2,
                 sagittal_size                      = 0.01,
                 grazing_angle                      = 0.003,
                 orientation_of_reflection_plane    = Orientation.UP,
                 invert_tangent_component           = False,
                 height_profile_data_file           = "mirror.dat",
                 height_profile_data_file_dimension = 1,
                 height_amplification_coefficient   = 1.0,
                 diffraction_order                  = 1,
                 grooving_density_0                 =800, # groove density [lines/mm] (coefficient a0 in the polynomial groove density: a0 + a1*y + a2*y^2 + a3*y^3 + a4*y^4)
                 grooving_density_1                 =0.0, # groove density polynomial coefficient a1 [lines/mm^2]
                 grooving_density_2                 =0.0, # groove density polynomial coefficient a2 [lines/mm^3]
                 grooving_density_3                 =0.0, # groove density polynomial coefficient a3 [lines/mm^4]
                 grooving_density_4                 =0.0, # groove density polynomial coefficient a4 [lines/mm^5]
                 grooving_angle                     = 0.0 # angle between the grove direction and the sagittal direction of the substrate
                 ):

        Grating.__init__(self,
                        name=name,
                        boundary_shape=Rectangle(x_left=-0.5*(sagittal_size),
                                                 x_right=0.5*(sagittal_size),
                                                 y_bottom=-0.5*tangential_size,
                                                 y_top=0.5*tangential_size),
                        surface_shape=self.get_shape(),
                        ruling=grooving_density_0*1e3)

        self.tangential_size                                  = tangential_size
        self.sagittal_size                                    = sagittal_size
        self.grazing_angle                                    = grazing_angle
        self.orientation_of_reflection_plane                  = orientation_of_reflection_plane
        self.invert_tangent_component                         = invert_tangent_component

        self.height_profile_data_file = height_profile_data_file
        self.height_profile_data_file_dimension = height_profile_data_file_dimension
        self.height_amplification_coefficient = height_amplification_coefficient

        self.diffraction_order = diffraction_order

        self.grooving_density_0 = grooving_density_0
        self.grooving_density_1 = grooving_density_1
        self.grooving_density_2 = grooving_density_2
        self.grooving_density_3 = grooving_density_3
        self.grooving_density_4 = grooving_density_4
        self.grooving_angle     = grooving_angle

    def get_alpha_angle(self):
        return self.grazing_angle

    def get_beta_angle(self, photon_energy):
        wavelength = 1.239842e-3/photon_energy # in mm

        return numpy.asin(wavelength*self.grooving_density_0 - numpy.cos(self.get_alpha_angle())) # Grating Output Angle

    def get_deflection_angle(self, photon_energy):
        return self.get_alpha_angle() + self.get_beta_angle(photon_energy) + 1.57079632679 # Grating Deflection Angle

    def get_output_orientation_vectors(self, photon_energy):
        deflection_angle = self.get_deflection_angle(photon_energy)
        tangent = 1 if not self.invert_tangent_component else -1

        if self.orientation_of_reflection_plane == Orientation.UP:
            return 0,  -numpy.sin(deflection_angle),  numpy.cos(deflection_angle),  0,  tangent
        elif self.orientation_of_reflection_plane == Orientation.DOWN:
            return 0,  numpy.sin(deflection_angle),  numpy.cos(deflection_angle),  0,  tangent
        elif self.orientation_of_reflection_plane == Orientation.LEFT:
            return numpy.sin(deflection_angle),  0, numpy.cos(deflection_angle),  tangent,  0
        elif self.orientation_of_reflection_plane == Orientation.RIGHT:
            return -numpy.sin(deflection_angle),  0, numpy.cos(deflection_angle),  tangent,  0

    def get_shape(self):
        raise NotImplementedError()

    def applyOpticalElement(self, wavefront, parameters=None):
        if not parameters.has_additional_parameter("srw_oe_wavefront_propagation_parameters"):
            wavefront_propagation_parameters = WavefrontPropagationParameters()
        else:
            wavefront_propagation_parameters = parameters.get_additional_parameter("srw_oe_wavefront_propagation_parameters")

            if not isinstance(wavefront_propagation_parameters, WavefrontPropagationParameters):
                raise ValueError("SRW Wavefront Propagation Parameters not present")

        substrate_mirror = self.get_substrate_mirror()

        grating = SRWLOptG(_mirSub=substrate_mirror,
                           _m=self.diffraction_order,
                           _grDen =self.grooving_density_0,
                           _grDen1=self.grooving_density_1,
                           _grDen2=self.grooving_density_2,
                           _grDen3=self.grooving_density_3,
                           _grDen4=self.grooving_density_4,
                           _grAng=self.grooving_angle)

        wavefront_propagation_parameters.to_SRW_array()

        optical_elements = [grating]
        propagation_parameters = [wavefront_propagation_parameters.to_SRW_array()]

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
                                                        _ang_r=self.get_deflection_angle(wavefront.get_photon_energy()),
                                                        _dim=dim,
                                                        _amp_coef=self.height_amplification_coefficient)

            elif self.height_profile_data_file_dimension == 2:
                height_profile_data = srwl_uti_read_data_cols(self.height_profile_data_file,
                                                              _str_sep='\t')

                optTrEr = srwl_opt_setup_surf_height_2d(_height_prof_data=height_profile_data,
                                                        _ang=self.grazing_angle,
                                                        _ang_r=self.get_deflection_angle(wavefront.get_photon_energy()),
                                                        _dim=dim,
                                                        _amp_coef=self.height_amplification_coefficient)

            optical_elements.append([optTrEr])
            propagation_parameters.append([WavefrontPropagationParameters().to_SRW_array()])


        optBL = SRWLOptC(optical_elements, propagation_parameters)

        srwl.PropagElecField(wavefront, optBL)

        return wavefront

    def get_substrate_mirror(self):
        nvx, nvy, nvz, tvx, tvy = self.get_orientation_vectors()
        x, y = self.getXY()

        return self.get_SRWLOptMir(nvx, nvy, nvz, tvx, tvy, x, y)

    def get_SRWLOptMir(self, nvx, nvy, nvz, tvx, tvy, x, y):
        mirror = SRWLOptMir()

        if isinstance(self.boundary_shape, Rectangle):
            ap_shape = ApertureShape.RECTANGULAR
        elif isinstance(self.boundary_shape, Ellipse) or  isinstance(self.boundary_shape, Circle):
            ap_shape = ApertureShape.ELLIPTIC

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
