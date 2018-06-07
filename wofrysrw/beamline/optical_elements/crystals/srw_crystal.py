import numpy

from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry
from syned.beamline.shape import Ellipse, Rectangle, Circle, Plane

from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement
from wofrysrw.propagator.wavefront2D.srw_wavefront import WavefrontPropagationParameters

from srwlib import SRWLOptCryst
from srwlib import srwl, srwl_opt_setup_surf_height_1d, srwl_opt_setup_surf_height_2d, srwl_uti_read_data_cols

from wofrysrw.beamline.optical_elements.mirrors.srw_mirror import Orientation, TreatInputOutput, ApertureShape, SimulationMethod

'''
        :param _d_sp: (_d_space) crystal reflecting planes d-spacing (John's dA) [A]
        :param _psi0r: real part of 0-th Fourier component of crystal polarizability (John's psi0c.real) (units?)
        :param _psi0i: imaginary part of 0-th Fourier component of crystal polarizability (John's psi0c.imag) (units?)
        :param _psi_hr: (_psiHr) real part of H-th Fourier component of crystal polarizability (John's psihc.real) (units?)
        :param _psi_hi: (_psiHi) imaginary part of H-th Fourier component of crystal polarizability (John's psihc.imag) (units?)
        :param _psi_hbr: (_psiHBr:) real part of -H-th Fourier component of crystal polarizability (John's psimhc.real) (units?)
        :param _psi_hbi: (_psiHBi:) imaginary part of -H-th Fourier component of crystal polarizability (John's psimhc.imag) (units?)
        :param _tc: crystal thickness [m] (John's thicum)
        :param _ang_as: (_Tasym) asymmetry angle [rad] (John's alphdg)
        :param _nvx: horizontal coordinate of outward normal to crystal surface (John's angles: thdg, chidg, phidg)
        :param _nvy: vertical coordinate of outward normal to crystal surface (John's angles: thdg, chidg, phidg)
        :param _nvz: longitudinal coordinate of outward normal to crystal surface (John's angles: thdg, chidg, phidg)
        :param _tvx: horizontal coordinate of central tangential vector (John's angles: thdg, chidg, phidg)
        :param _tvy: vertical coordinate of central tangential vector (John's angles: thdg, chidg, phidg)
        :param _uc: crystal use case: 1- Bragg Reflection, 2- Bragg Transmission (Laue cases to be added)

'''

class SRWCrystal(Crystal, SRWOpticalElement):
    def __init__(self,
                 name                                 = "Undefined",
                 tangential_size                      = 1.2,
                 sagittal_size                        = 0.01,
                 vertical_position_of_mirror_center   = 0.0,
                 horizontal_position_of_mirror_center = 0.0,
                 orientation_of_reflection_plane      = Orientation.UP,
                 invert_tangent_component             = False,
                 d_spacing                            = 0.0,
                 psi_0r                               = 0.0,
                 psi_0i                               = 0.0,
                 psi_hr                               = 0.0,
                 psi_hi                               = 0.0,
                 psi_hbr                              = 0.0,
                 psi_hbi                              = 0.0,
                 asymmetry_angle                      = 0.0,
                 thickness                            = 0.0,
                 diffraction_geometry                 = DiffractionGeometry.BRAGG
                ):
        Crystal.__init__(name,
                         surface_shape=Plane(),
                         boundary_shape=Rectangle(x_left=horizontal_position_of_mirror_center - 0.5*sagittal_size,
                                                  x_right=horizontal_position_of_mirror_center + 0.5*sagittal_size,
                                                  y_bottom=vertical_position_of_mirror_center - 0.5*tangential_size,
                                                  y_top=vertical_position_of_mirror_center + 0.5*tangential_size),
                         material="Unknown",
                         diffraction_geometry=diffraction_geometry,
                         asymmetry_angle = asymmetry_angle,
                         thickness = thickness
                        )
        self.tangential_size                                  = tangential_size
        self.sagittal_size                                    = sagittal_size
        self.orientation_of_reflection_plane                  = orientation_of_reflection_plane
        self.invert_tangent_component                         = invert_tangent_component

        self.d_spacing                            = d_spacing
        self.psi_0r                               = psi_0r
        self.psi_0i                               = psi_0i
        self.psi_hr                               = psi_hr
        self.psi_hi                               = psi_hi
        self.psi_hbr                              = psi_hbr
        self.psi_hbi                              = psi_hbi
        self.asymmetry_angle                      = asymmetry_angle
        self.thickness                            = thickness
        self.diffraction_geometry                 = diffraction_geometry

        if diffraction_geometry == DiffractionGeometry.LAUE: raise NotImplementedError("Laue Geometry is not yet supported")

    def toSRWLOpt(self):
        nvx, nvy, nvz, tvx, tvy = self.get_orientation_vectors()

        return SRWLOptCryst(_d_sp=self.d_spacing,
                            _psi0r=self.psi_0r,
                            _psi0i=self.psi_0i,
                            _psi_hr=self.psi_hr,
                            _psi_hi=self.psi_hi,
                            _psi_hbr=self.psi_hbr,
                            _psi_hbi=self.psi_hbi,
                            _tc=self.thickness,
                            _ang_as=self.asymmetry_angle,
                            _nvx=nvx,
                            _nvy=nvy,
                            _nvz=nvz,
                            _tvx=tvx,
                            _tvy=tvy,
                            _uc=1 if self.diffraction_geometry==DiffractionGeometry.BRAGG else 0)

