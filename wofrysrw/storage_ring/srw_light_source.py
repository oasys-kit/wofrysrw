import numpy
import array

from srwlib import srwl, SRWLRadMesh, SRWLStokes

from syned.storage_ring.light_source import LightSource

from wofry.elements.decorators import WOLightSourceDecorator
from wofrysrw.storage_ring.srw_magnetic_structure import SRWMagneticStructure
from wofrysrw.storage_ring.srw_electron_beam import SRWElectronBeam, SRWElectronBeamGeometricalProperties
from wofrysrw.propagator.wavefront2D.srw_wavefront import SRWWavefront

class WavefrontPrecisionParameters(object):
    def __init__(self,
                 sr_method = 1,  #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
                 relative_precision = 0.01, # relative precision
                 start_integration_longitudinal_position = 0, # longitudinal position to start integration (effective if < zEndInteg)
                 end_integration_longitudinal_position = 0, # longitudinal position to finish integration (effective if > zStartInteg)
                 number_of_points_for_trajectory_calculation = 50000, #Number of points for trajectory calculation
                 use_terminating_terms = 1, # Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
                 sampling_factor_for_adjusting_nx_ny = 0.0 # sampling factor for adjusting nx, ny (effective if > 0)
                 ):
        self._sr_method = sr_method
        self._relative_precision = relative_precision
        self._start_integration_longitudinal_position = start_integration_longitudinal_position
        self._end_integration_longitudinal_position = end_integration_longitudinal_position
        self._number_of_points_for_trajectory_calculation = number_of_points_for_trajectory_calculation
        self._use_terminating_terms = use_terminating_terms
        self._sampling_factor_for_adjusting_nx_ny = sampling_factor_for_adjusting_nx_ny

    def to_SRW_array(self):
        return [int(self._sr_method),
                float(self._relative_precision),
                float(self._start_integration_longitudinal_position),
                float(self._end_integration_longitudinal_position),
                int(self._number_of_points_for_trajectory_calculation),
                int(self._use_terminating_terms),
                float(self._sampling_factor_for_adjusting_nx_ny)]


class PowerDensityPrecisionParameters(object):
    def __init__(self,
                 precision_factor = 1.5,
                 computation_method = 1, # (1- "near field", 2- "far field")
                 initial_longitudinal_position = 0.0, # (effective if initial_longitudinal_position < final_longitudinal_position)
                 final_longitudinal_position = 0.0, # (effective if initial_longitudinal_position < final_longitudinal_position)
                 number_of_points_for_trajectory_calculation = 20000 #number of points for (intermediate) trajectory calculation
                 ):
        self._precision_factor = precision_factor
        self._computation_method = computation_method
        self._initial_longitudinal_position = initial_longitudinal_position
        self._final_longitudinal_position = final_longitudinal_position
        self._number_of_points_for_trajectory_calculation = number_of_points_for_trajectory_calculation

    def to_SRW_array(self):
        return [float(self._precision_factor),
                int(self._computation_method),
                float(self._initial_longitudinal_position),
                float(self._final_longitudinal_position),
                int(self._number_of_points_for_trajectory_calculation)]

class SourceWavefrontParameters(object):
    def __init__(self, 
                 photon_energy_min = 100,
                 photon_energy_max = 10100,
                 photon_energy_points = 51,
                 h_slit_gap = 0,
                 h_slit_points = 1,
                 v_slit_gap = 0,
                 v_slit_points = 1,
                 distance = 10.0,
                 wavefront_precision_parameters=WavefrontPrecisionParameters()):
        self._photon_energy_min         = photon_energy_min
        self._photon_energy_max         = photon_energy_max
        self._photon_energy_points      = photon_energy_points
        self._h_slit_gap                = h_slit_gap
        self._h_slit_points             = h_slit_points
        self._v_slit_gap                = v_slit_gap
        self._v_slit_points             = v_slit_points
        self._distance                  = distance
        self._wavefront_precision_parameters  = wavefront_precision_parameters

    def to_SRWRadMesh(self):
        return SRWLRadMesh(self._photon_energy_min,
                           self._photon_energy_max,
                           self._photon_energy_points,
                           -self._h_slit_gap/2, self._h_slit_gap/2, self._h_slit_points,
                           -self._v_slit_gap/2, self._v_slit_gap/2, self._v_slit_points,
                           self._distance)

    def to_SRWLStokes(self):
        stk = SRWLStokes()
        stk.allocate(self._photon_energy_points,
                     self._h_slit_points,
                     self._v_slit_points)
        stk.mesh = self.to_SRWRadMesh()

        return stk

class PolarizationComponent:
    LINEAR_HORIZONTAL  = 0
    LINEAR_VERTICAL    = 1
    LINEAR_45_DEGREES  = 2
    LINEAR_135_DEGREES = 3
    CIRCULAR_RIGHT     = 4
    CIRCULAR_LEFT      = 5
    TOTAL              = 6

'''
:param polarization_component_to_be_extracted:
               =0 -Linear Horizontal;
               =1 -Linear Vertical;
               =2 -Linear 45 degrees;
               =3 -Linear 135 degrees;
               =4 -Circular Right;
               =5 -Circular Left;
               =6 -Total
:param calculation_type:
               =0 -"Single-Electron" Intensity;
               =1 -"Multi-Electron" Intensity;
               =2 -"Single-Electron" Flux;
               =3 -"Multi-Electron" Flux;
               =4 -"Single-Electron" Radiation Phase;
               =5 -Re(E): Real part of Single-Electron Electric Field;
               =6 -Im(E): Imaginary part of Single-Electron Electric Field;
               =7 -"Single-Electron" Intensity, integrated over Time or Photon Energy (i.e. Fluence)
:param type_of_dependence: 
               =0 -vs e (photon energy or time);
               =1 -vs x (horizontal position or angle);
               =2 -vs y (vertical position or angle);
               =3 -vs x&y (horizontal and vertical positions or angles);
               =4 -vs e&x (photon energy or time and horizontal position or angle);
               =5 -vs e&y (photon energy or time and vertical position or angle);
               =6 -vs e&x&y (photon energy or time, horizontal and vertical positions or angles);
:param fixed_input_photon_energy_or_time: input photon energy [eV] or time [s] to keep fixed (to be taken into account for dependences vs x, y, x&y)
:param fixed_horizontal_position: input horizontal position [m] to keep fixed (to be taken into account for dependences vs e, y, e&y)
:param fixed_vertical_position: input vertical position [m] to keep fixed (to be taken into account for dependences vs e, x, e&x)
'''
class FluxCalculationParameters(object):
    def __init__(self,
                 polarization_component_to_be_extracted=6,
                 calculation_type=0,
                 type_of_dependence=0,
                 fixed_input_photon_energy_or_time = 0.0,
                 fixed_horizontal_position = 0.0,
                 fixed_vertical_position = 0.0):
    
        self._polarization_component_to_be_extracted = polarization_component_to_be_extracted
        self._calculation_type                       = calculation_type                      
        self._type_of_dependence                     = type_of_dependence                    
        self._fixed_input_photon_energy_or_time      = fixed_input_photon_energy_or_time     
        self._fixed_horizontal_position              = fixed_horizontal_position             
        self._fixed_vertical_position                = fixed_vertical_position                

class PhotonSourceProperties(object):

    def __init__(self,
                 rms_h = 0.0,
                 rms_v = 0.0,
                 rms_hp = 0.0,
                 rms_vp = 0.0,
                 coherence_volume_h = 0.0,
                 coherence_volume_v = 0.0,
                 diffraction_limit = 0.0):
        self._rms_h = rms_h
        self._rms_v = rms_v
        self._rms_hp = rms_hp
        self._rms_vp = rms_vp
        self._coherence_volume_h = coherence_volume_h
        self._coherence_volume_v = coherence_volume_v
        self._diffraction_limit = diffraction_limit

    def to_info(self):
        info = 'Photon beam (convolution): \n'
        info += '   RMS size H/V [um]: '+ repr(self._rms_h*1e6) + '  /  ' + repr(self._rms_v*1e6) + '\n'
        info += '   RMS divergence H/V [urad]: '+ repr(self._rms_hp*1e6) + '  /  ' + repr(self._rms_vp*1e6) + '\n\n'
        info += '   Coherent volume in H phase space: '+ repr(self._coherence_volume_h) + '\n'
        info += '   Coherent volume in V phase space: '+ repr(self._coherence_volume_v) + '\n\n'
        info += '   RMS diffraction limit source size [um]: '+ repr(self._diffraction_limit*1e6) + '\n'
        info += '   FWHM diffraction limit source size [um]: '+ repr(self._diffraction_limit*2.35*1e6)

        return info

class SRWLightSource(LightSource, WOLightSourceDecorator):
    def __init__(self,
                 name="Undefined",
                 electron_energy_in_GeV = 1.0,
                 electron_energy_spread = 0.0,
                 ring_current = 0.1,
                 number_of_bunches = 400,
                 electron_beam_size_h=1e-5,
                 electron_beam_size_v=1e-5,
                 electron_beam_divergence_h=0.0,
                 electron_beam_divergence_v=0.0,
                 magnetic_structure=SRWMagneticStructure()):
        electron_beam = SRWElectronBeam(energy_in_GeV=electron_energy_in_GeV,
                                        energy_spread=electron_energy_spread,
                                        current=ring_current,
                                        number_of_bunches=number_of_bunches)

        electron_beam.set_moments_from_electron_beam_geometrical_properties(SRWElectronBeamGeometricalProperties(electron_beam_size_h=electron_beam_size_h,
                                                                                                                 electron_beam_divergence_h=electron_beam_divergence_h,
                                                                                                                 electron_beam_size_v=electron_beam_size_v,
                                                                                                                 electron_beam_divergence_v=electron_beam_divergence_v))

        LightSource.__init__(self, name, electron_beam, magnetic_structure)

    def get_gamma(self):
        return self._electron_beam.gamma()

    def get_photon_source_properties(self):
        return NotImplementedError("must be implemented in subclasses")

    # from Wofry Decorator
    def get_wavefront(self, wavefront_parameters):
        return self.get_SRW_Wavefront(source_wavefront_parameters=wavefront_parameters).toGenericWavefront()

    def get_SRW_Wavefront(self, source_wavefront_parameters = SourceWavefrontParameters()):
        mesh = source_wavefront_parameters.to_SRWRadMesh()

        wfr = SRWWavefront()
        wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
        wfr.mesh = mesh
        wfr.partBeam = self._electron_beam.to_SRWLPartBeam()

        srwl.CalcElecFieldSR(wfr,
                             0,
                             self._magnetic_structure.get_SRWLMagFldC(),
                             source_wavefront_parameters._wavefront_precision_parameters.to_SRW_array())

        return wfr

    def get_intensity_from_electric_field(self,
                                          output_array,
                                          srw_wavefront,
                                          flux_calculation_parameters = FluxCalculationParameters()):
        srwl.CalcIntFromElecField(output_array,
                                  srw_wavefront,
                                  flux_calculation_parameters._polarization_component_to_be_extracted,
                                  flux_calculation_parameters._calculation_type,
                                  flux_calculation_parameters._type_of_dependence,
                                  flux_calculation_parameters._fixed_input_photon_energy_or_time,
                                  flux_calculation_parameters._fixed_horizontal_position,
                                  flux_calculation_parameters._fixed_vertical_position)

        return output_array


    def get_flux_per_unit_surface(self, source_wavefront_parameters = SourceWavefrontParameters(), multi_electron=True):

        flux_calculation_parameters=FluxCalculationParameters(calculation_type                  = 1 if multi_electron==True else 0,
                                                              type_of_dependence                = 3)

        srw_wavefront = self.get_SRW_Wavefront(source_wavefront_parameters)

        h_array=numpy.linspace(srw_wavefront.mesh.xStart, srw_wavefront.mesh.xFin, srw_wavefront.mesh.nx)
        v_array=numpy.linspace(srw_wavefront.mesh.yStart, srw_wavefront.mesh.yFin, srw_wavefront.mesh.ny)
        e_array=numpy.linspace(srw_wavefront.mesh.eStart, srw_wavefront.mesh.eFin, srw_wavefront.mesh.ne)

        flux_per_unit_surface_array = numpy.zeros((e_array.size, h_array.size, v_array.size))

        for ie in range(e_array.size):
            output_array = array.array('f', [0]*srw_wavefront.mesh.nx*srw_wavefront.mesh.ny) #"flat" array to take 2D intensity data

            flux_calculation_parameters._fixed_input_photon_energy_or_time = e_array[ie]
            self.get_intensity_from_electric_field(output_array, srw_wavefront, flux_calculation_parameters)

            data = numpy.ndarray(buffer=output_array, shape=(srw_wavefront.mesh.ny, srw_wavefront.mesh.nx),dtype=output_array.typecode)

            for ix in range(h_array.size):
                for iy in range(v_array.size):
                    flux_per_unit_surface_array[ie,ix,iy,] = data[iy,ix]
        
        return (e_array, h_array, v_array, flux_per_unit_surface_array)

    def get_spectral_flux(self, source_wavefront_parameters = SourceWavefrontParameters(), multi_electron=True):

        flux_calculation_parameters=FluxCalculationParameters(calculation_type   = 1 if multi_electron else 0,
                                                              type_of_dependence = 0)

        srw_wavefront = self.get_SRW_Wavefront(source_wavefront_parameters)

        output_array = array.array('f', [0]*srw_wavefront.mesh.ne)

        self.get_intensity_from_electric_field(output_array, srw_wavefront, flux_calculation_parameters)
        
        data = numpy.ndarray(buffer=output_array, shape=srw_wavefront.mesh.ne, dtype=output_array.typecode)

        energy_array=numpy.linspace(srw_wavefront.mesh.eStart,
                                    srw_wavefront.mesh.eFin,
                                    srw_wavefront.mesh.ne)
        spectral_flux_array = numpy.zeros(energy_array.size)

        for ie in range(energy_array.size):
            spectral_flux_array[ie] = data[ie]

        return (energy_array, spectral_flux_array)



    def get_power_density(self,
                          source_wavefront_parameters = SourceWavefrontParameters(),
                          power_density_precision_parameters = PowerDensityPrecisionParameters()):

        stkP = source_wavefront_parameters.to_SRWLStokes()

        srwl.CalcPowDenSR(stkP,
                          self._electron_beam.to_SRWLPartBeam(),
                          0,
                          self._magnetic_structure.get_SRWLMagFldC(),
                          power_density_precision_parameters.to_SRW_array())

        hArray = numpy.zeros(stkP.mesh.nx)
        vArray = numpy.zeros(stkP.mesh.ny)
        powerArray = numpy.zeros((stkP.mesh.nx,stkP.mesh.ny))

        # fill arrays
        ij = -1
        for j in range(stkP.mesh.ny):
            for i in range(stkP.mesh.nx):
                ij += 1
                xx = stkP.mesh.xStart + i*(stkP.mesh.xFin-stkP.mesh.xStart)/(stkP.mesh.nx-1)
                yy = stkP.mesh.yStart + j*(stkP.mesh.yFin-stkP.mesh.yStart)/(stkP.mesh.ny-1)
                powerArray[i,j] = stkP.arS[ij]
                hArray[i] = xx # mm
                vArray[j] = yy # mm

        return (hArray, vArray, powerArray)

    @classmethod
    def get_total_power_from_power_density(cls, h_array, v_array, power_density_matrix):
        area = (numpy.abs(h_array[1]-h_array[0])*numpy.abs(v_array[1]-v_array[0]))*1e6
        total_power = 0
        for i in range(0, len(h_array)):
            for j in range(0, len(v_array)):
                total_power += power_density_matrix[i, j]*area

        return total_power
