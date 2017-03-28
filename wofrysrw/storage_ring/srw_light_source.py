import numpy
import array

import scipy.constants as codata

from srwlib import srwl, SRWLRadMesh

from syned.storage_ring.light_source import LightSource

from wofrysrw.storage_ring.srw_magnetic_structure import SRWMagneticStructure
from wofrysrw.storage_ring.srw_electron_beam import SRWElectronBeam, SRWElectronBeamGeometricalProperties
from wofrysrw.propagator.wavefront2D.srw_wavefront import SRWWavefront

class SRWPrecisionParameters(object):
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

class SourceWavefrontParameters(object):
    def __init__(self, 
                 photon_energy_min = 100,
                 photon_energy_max = 10100,
                 photon_energy_points = 51,
                 h_slit_gap = 1e-3, 
                 h_slit_points = 51, 
                 v_slit_gap = 1e-3, 
                 v_slit_points = 51, 
                 distance = 10,
                 srw_precision_parameters=SRWPrecisionParameters()):
        self._photon_energy_min         = photon_energy_min
        self._photon_energy_max         = photon_energy_max
        self._photon_energy_points      = photon_energy_points
        self._h_slit_gap                = h_slit_gap
        self._h_slit_points             = h_slit_points
        self._v_slit_gap                = v_slit_gap
        self._v_slit_points             = v_slit_points
        self._distance                  = distance
        self._srw_precision_parameters  = srw_precision_parameters

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

class SRWLightSource(LightSource):
    def __init__(self,
                 name="Undefined",
                 electron_energy_in_GeV = 1.0,
                 electron_energy_spread = 0.0,
                 ring_current = 0.1,
                 electrons_per_bunch = 400,
                 electron_beam_size_h=0.0,
                 electron_beam_size_v=0.0,
                 emittance=0.0,
                 coupling_costant=0.0,
                 magnetic_structure=SRWMagneticStructure()):
        electron_beam = SRWElectronBeam(energy_in_GeV=electron_energy_in_GeV,
                                          energy_spread=electron_energy_spread,
                                          current=ring_current,
                                          electrons_per_bunch=electrons_per_bunch)

        electron_beam.set_moments_from_electron_beam_geometrical_properties(SRWElectronBeamGeometricalProperties(electron_beam_size_h=electron_beam_size_h,
                                                                                                                   electron_beam_divergence_h=(emittance/electron_beam_size_h),
                                                                                                                   electron_beam_size_v=electron_beam_size_v,
                                                                                                                   electron_beam_divergence_v=(coupling_costant*emittance/electron_beam_size_v)))

        LightSource.__init__(self, name, electron_beam, magnetic_structure)

    def get_gamma(self):
        return self._electron_beam.gamma()

    def get_photon_source_properties(self):
        return NotImplementedError("must be implemented in subclasses")

    def get_SRW_Wavefront(self, source_wavefront_parameters = SourceWavefrontParameters()):

        mesh = SRWLRadMesh(source_wavefront_parameters._photon_energy_min,
                           source_wavefront_parameters._photon_energy_max,
                           source_wavefront_parameters._photon_energy_points,
                           -source_wavefront_parameters._h_slit_gap/2, source_wavefront_parameters._h_slit_gap/2, source_wavefront_parameters._h_slit_points,
                           -source_wavefront_parameters._v_slit_gap/2, source_wavefront_parameters._v_slit_gap/2, source_wavefront_parameters._v_slit_points,
                           source_wavefront_parameters._distance)

        wfr = SRWWavefront()
        wfr.mesh = mesh
        wfr.partBeam = self._electron_beam.to_SRWLPartBeam()
        wfr.allocate(mesh.ne, mesh.nx, mesh.ny)

        srwl.CalcElecFieldSR(wfr, 0, self._magnetic_structure.get_SRWLMagFldC(), source_wavefront_parameters._srw_precision_parameters.to_SRW_array())

        return wfr

    def get_radiation(self, srw_wavefront):

        mesh0 = srw_wavefront.mesh

        INTENSITY_TYPE_MULTI_ELECTRON=1

        hArray=numpy.linspace(srw_wavefront.mesh.xStart, srw_wavefront.mesh.xFin, srw_wavefront.mesh.nx)
        vArray=numpy.linspace(srw_wavefront.mesh.yStart, srw_wavefront.mesh.yFin, srw_wavefront.mesh.ny)
        eArray=numpy.linspace(srw_wavefront.mesh.eStart, srw_wavefront.mesh.eFin, srw_wavefront.mesh.ne)

        intensArray = numpy.zeros((eArray.size,hArray.size,vArray.size,))
        for ie in range(eArray.size):
            arI0 = array.array('f', [0]*mesh0.nx*mesh0.ny) #"flat" array to take 2D intensity data
            srwl.CalcIntFromElecField(arI0, srw_wavefront, 6, INTENSITY_TYPE_MULTI_ELECTRON, 3, eArray[ie], 0, 0) # 6 is for total polarizarion; 0=H, 1=V

            data = numpy.ndarray(buffer=arI0, shape=(mesh0.ny, mesh0.nx),dtype=arI0.typecode)

            for ix in range(hArray.size):
                for iy in range(vArray.size):
                    intensArray[ie,ix,iy,] = data[iy,ix]
        
        return (eArray, hArray, vArray, intensArray)

    @classmethod
    def get_spectrum_from_radiation(cls, h_array, v_array, radiation_matrix):
        return (radiation_matrix.sum(axis=2)).sum(axis=1)*(h_array[1]-h_array[0])*(v_array[1]-v_array[0])

    @classmethod
    def get_power_density_from_radiation(cls, energy_array, radiation_matrix):
        return radiation_matrix.sum(axis=0)*(energy_array[1]-energy_array[0])*codata.e*1e3

    @classmethod
    def get_total_power_from_power_density(cls, h_array, v_array, power_density_matrix):
        area = numpy.abs(h_array[1]-h_array[0])*numpy.abs(v_array[1]-v_array[0])
        total_power = 0
        for i in range(0, len(h_array)):
            for j in range(0, len(v_array)):
                total_power += power_density_matrix[i, j]*area

        return total_power
