import numpy as np

from syned.beamline.optical_elements.crystals.crystal import Crystal, DiffractionGeometry
from syned.beamline.shape import  Rectangle, Plane

from wofrysrw.beamline.optical_elements.srw_optical_element import SRWOpticalElement
from wofrysrw.beamline.optical_elements.mirrors.srw_mirror import Orientation
from wofrysrw.util.srw import SRWLOptCryst, pi, srwl_uti_cryst_pol_f, srwl_uti_cryst_pl_sp


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
                 name: str                                 = "Undefined",
                 optical_element_displacement  = None,
                 orientation_of_reflection_plane: int      = Orientation.UP,
                 material                                  = None,
                 miller_indices                    = [1, 1, 1],
                 d_spacing:float                      = 0.0,
                 psi_0r:float                         = 0.0,
                 psi_0i:float                         = 0.0,
                 psi_hr:float                         = 0.0,
                 psi_hi:float                         = 0.0,
                 psi_hbr:float                        = 0.0,
                 psi_hbi:float                        = 0.0,
                 asymmetry_angle:float                = 0.0,
                 thickness:float                      = 0.0,
                 diffraction_geometry:int                  = DiffractionGeometry.BRAGG,
                 energy:float                         = 0.0
                ):
        SRWOpticalElement.__init__(self, optical_element_displacement=optical_element_displacement)

        Crystal.__init__(self,
                         name,
                         surface_shape=Plane(),
                         boundary_shape=Rectangle(x_left=0.0,
                                                  x_right=0.0,
                                                  y_bottom=0.0,
                                                  y_top=0.0),
                         material="Unknown",
                         diffraction_geometry=diffraction_geometry,
                         asymmetry_angle = asymmetry_angle,
                         thickness = thickness
                        )

        self.orientation_of_reflection_plane = orientation_of_reflection_plane

        self.material       = material
        self.miller_indices = miller_indices

        self.d_spacing            = d_spacing
        self.psi_0r               = psi_0r
        self.psi_0i               = psi_0i
        self.psi_hr               = psi_hr
        self.psi_hi               = psi_hi
        self.psi_hbr              = psi_hbr
        self.psi_hbi              = psi_hbi
        self.asymmetry_angle      = asymmetry_angle
        self.thickness            = thickness
        self.diffraction_geometry = diffraction_geometry
        self.energy               = energy

        if diffraction_geometry == DiffractionGeometry.LAUE: raise NotImplementedError("Laue Geometry is not yet supported")

        if not material is None:
            self.d_spacing                 = srwl_uti_cryst_pl_sp(_hr=self.miller_indices, _mat=self.material)
            psi_0r, psi_0i, psi_hr, psi_hi = srwl_uti_cryst_pol_f(_en=self.energy,_hr=self.miller_indices,_mat=self.material)

            self.psi_0r  = psi_0r
            self.psi_0i  = psi_0i
            self.psi_hr  = psi_hr
            self.psi_hi  = psi_hi
            self.psi_hbr = psi_hr
            self.psi_hbi = psi_hi

        self._srw_object =  SRWLOptCryst(_d_sp=self.d_spacing,
                                         _psi0r=self.psi_0r,
                                         _psi0i=self.psi_0i,
                                         _psi_hr=self.psi_hr,
                                         _psi_hi=self.psi_hi,
                                         _psi_hbr=self.psi_hbr,
                                         _psi_hbi=self.psi_hbi,
                                         _tc=self.thickness,
                                         _ang_as=self.asymmetry_angle,
                                         _uc=1 if self.diffraction_geometry == DiffractionGeometry.BRAGG else 0)

        if   orientation_of_reflection_plane == Orientation.LEFT:  orientation_angle = 0.5*pi
        elif orientation_of_reflection_plane == Orientation.RIGHT: orientation_angle = 1.5*pi
        elif orientation_of_reflection_plane == Orientation.UP:    orientation_angle = 0.0
        elif orientation_of_reflection_plane == Orientation.DOWN:  orientation_angle = pi

        self._orientation_data = self._srw_object.find_orient(self.energy, orientation_angle)

        nvx, nvy, nvz, tvx, tvy = self.get_orientation_vectors()

        self._srw_object.set_orient(nvx, nvy, nvz, tvx, tvy)

        if   orientation_of_reflection_plane == Orientation.LEFT:  self.grazing_angle = -np.arccos(nvx) + np.pi
        elif orientation_of_reflection_plane == Orientation.RIGHT: self.grazing_angle = np.arccos(nvx)
        elif orientation_of_reflection_plane == Orientation.UP:    self.grazing_angle = np.arccos(nvy)
        elif orientation_of_reflection_plane == Orientation.DOWN:  self.grazing_angle = -np.arccos(nvy) + np.pi

    def get_orientation_vectors(self):
        vectors_data = self._orientation_data[0]
        tO = vectors_data[0]
        nO = vectors_data[2]

        return nO[0], nO[1], nO[2], tO[0], tO[1]

    def get_output_orientation_vectors(self):
        output_orientation = self._orientation_data[1]
        rxO = output_orientation[0]
        rzO = output_orientation[2]

        return rzO[0], rzO[1], rzO[2], rxO[0], rxO[1]

    def toSRWLOpt(self):
        return self._srw_object

    def to_python_code(self, data=None):
        oe_name = data[0]

        text_code  = oe_name + "="+ "SRWLOptCryst(_d_sp=" + str(self.d_spacing) + "," + "\n"
        text_code += "                        _psi0r=" + str(self.psi_0r) + "," + "\n"
        text_code += "                        _psi0i=" + str(self.psi_0i) + "," + "\n"
        text_code += "                        _psi_hr=" + str(self.psi_hr) + "," + "\n"
        text_code += "                        _psi_hi=" + str(self.psi_hi) + "," + "\n"
        text_code += "                        _psi_hbr=" + str(self.psi_hbr) + "," + "\n"
        text_code += "                        _psi_hbi=" + str(self.psi_hbi) + "," + "\n"
        text_code += "                        _tc=" + str(self.thickness) + "," + "\n"
        text_code += "                        _ang_as=" + str(self.asymmetry_angle) + "," + "\n"
        text_code += "                        _uc=" + str(1 if self.diffraction_geometry==DiffractionGeometry.BRAGG else 0) + ")" + "\n"

        text_code += f"\norientation_data = {oe_name}.find_orient({self.energy}, 0.)\n"
        text_code += f"tO = orientation_data[0][0]\n"
        text_code += f"nO = orientation_data[0][2]\n"
        text_code += f"{oe_name}.set_orient(nO[0], nO[1], nO[2], tO[0], tO[1])\n"

        return text_code
