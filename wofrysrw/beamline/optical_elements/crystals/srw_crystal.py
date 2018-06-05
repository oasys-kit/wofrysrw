import numpy

from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry
from syned.beamline.shape import Ellipse, Rectangle, Circle

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
                 name                               = "Undefined",
                 diffraction_geometry=DiffractionGeometry.BRAGG,
                 miller_index_h=1,
                 miller_index_k=1,
                 miller_index_l=1,
                 asymmetry_angle=0.0,
                 thickness=0.0,
                 psi0r=0.0,
                 psi0i=0.0,
                 psi_hr=0.0,
                 psi_hi=0.0,
                 psi_hbr=0.0,
                 psi_hbi=0.0):
        pass
