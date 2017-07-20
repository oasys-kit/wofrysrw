import numpy

from srwlib import srwl

from syned.storage_ring.light_source import LightSource

from wofry.beamline.decorators import WOLightSourceDecorator
from wofrysrw.storage_ring.srw_magnetic_structure import SRWMagneticStructure
from wofrysrw.storage_ring.srw_electron_beam import SRWElectronBeam, SRWElectronBeamGeometricalProperties
from wofrysrw.propagator.wavefront2D.srw_wavefront import SRWWavefront, WavefrontParameters

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

    def get_SRW_Wavefront(self, source_wavefront_parameters = WavefrontParameters()):
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

    def get_intensity(self, source_wavefront_parameters = WavefrontParameters(), multi_electron=True):
        
        srw_wavefront = self.get_SRW_Wavefront(source_wavefront_parameters)
        
        return srw_wavefront.get_intensity(multi_electron=multi_electron)

    def get_flux(self, source_wavefront_parameters = WavefrontParameters(), multi_electron=True):

        srw_wavefront = self.get_SRW_Wavefront(source_wavefront_parameters)

        return srw_wavefront.get_flux(multi_electron=multi_electron)

    def get_power_density(self,
                          source_wavefront_parameters = WavefrontParameters(),
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