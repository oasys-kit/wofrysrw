from wofrysrw.beamline.optical_elements.mirrors.srw_mirror import SRWMirror, Orientation, SimulationMethod, TreatInputOutput
from syned.beamline.shape import ParabolicCylinder, Side, Convexity, Paraboloid

from wofrysrw.util.srw import SRWLOptMirPar

class SRWParaboloidMirror(SRWMirror):
    def __init__(self,
                 name :str                                  = "Undefined",
                 optical_element_displacement               = None,
                 tangential_size :float                = 1.2,
                 sagittal_size :float                  = 0.01,
                 grazing_angle:float                   = 0.003,
                 orientation_of_reflection_plane :int       = Orientation.UP,
                 invert_tangent_component :bool             = False,
                 focal_length :float                   = 0.0,
                 at_infinity : int                          = Side.SOURCE,
                 sagittal_radius                      = 1.e23,
                 height_profile_data_file :str              = "mirror.dat",
                 height_profile_data_file_dimension :int     = 1,
                 height_amplification_coefficient :float= 1.0):
        
        if sagittal_radius > 1e8: shape = ParabolicCylinder(parabola_parameter=2*focal_length, at_infinity=at_infinity, convexity=Convexity.UPWARD)
        else:                     shape = Paraboloid(parabola_parameter=2*focal_length, at_infinity=at_infinity)
        
        self._sagittal_radius = 1.e23 if sagittal_radius > 1e8 else sagittal_radius
        
        super().__init__(name=name,
                         shape=shape,
                         optical_element_displacement=optical_element_displacement,
                         tangential_size=tangential_size,
                         sagittal_size=sagittal_size,
                         grazing_angle=grazing_angle,
                         orientation_of_reflection_plane=orientation_of_reflection_plane,
                         invert_tangent_component=invert_tangent_component,
                         height_profile_data_file=height_profile_data_file,
                         height_profile_data_file_dimension=height_profile_data_file_dimension,
                         height_amplification_coefficient=height_amplification_coefficient)
    
    def get_SRWLOptMir(self, nvx, nvy, nvz, tvx, tvy, x, y, ap_shape):
        return SRWLOptMirPar(_f=0.5*self.get_surface_shape().get_parabola_parameter(),
                             _uc='f' if self.get_surface_shape().get_at_infinity() == Side.SOURCE else 'c',
                             _r_sag=self._sagittal_radius,
                             _size_tang=self.tangential_size,
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

    def fromSRWLOpt(self, srwlopt: SRWLOptMirPar):
        if not isinstance(srwlopt, SRWLOptMirPar):
            raise ValueError("SRW object is not a SRWLOptMirPar object")
        
        at_infinity        = Side.SOURCE if srwlopt.uc == 'f' else Side.IMAGE
        parabola_parameter = 2*srwlopt.f
        
        if srwlopt.radSag > 1e8: shape = ParabolicCylinder(parabola_parameter=parabola_parameter, at_infinity=at_infinity, convexity=Convexity.UPWARD)
        else:                    shape = Paraboloid(parabola_parameter=parabola_parameter, at_infinity=at_infinity)

        super().fromSRWLOpt(srwlopt, shape)

        self._sagittal_radius = 1.e23 if srwlopt.radSag > 1e8 else srwlopt.radSag

    def to_python_code_aux(self, nvx, nvy, nvz, tvx, tvy, x, y, ap_shape):
        text_code  = f"SRWLOptMirPar(_f={0.5*self.get_surface_shape().get_parabola_parameter()},\n"
        text_code += f"              _uc='{'f' if self.get_surface_shape().get_at_infinity() == Side.SOURCE else 'c'}',\n"
        text_code += f"              _r_sag={self._sagittal_radius},\n"
        text_code += f"              _size_tang={self.tangential_size},\n"
        text_code += f"              _size_sag={self.sagittal_size},\n"
        text_code += f"              _ap_shape='{ap_shape}',\n"
        text_code += f"              _sim_meth={SimulationMethod.THICK},\n"
        text_code += f"              _treat_in_out={TreatInputOutput.WAVEFRONT_INPUT_CENTER_OUTPUT_CENTER},\n"
        text_code += f"              _nvx={nvx},\n"
        text_code += f"              _nvy={nvy},\n"
        text_code += f"              _nvz={nvz},\n"
        text_code += f"              _tvx={tvx},\n"
        text_code += f"              _tvy={tvy},\n"
        text_code += f"              _x={x},\n"
        text_code += f"              _y={y})\n"

        return text_code
