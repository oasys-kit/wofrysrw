import numpy

from srwlib import srwl

from wofrysrw.storage_ring.srw_light_source import SRWLightSource, PhotonSourceProperties, SourceWavefrontParameters
from wofrysrw.storage_ring.magnetic_structures.srw_undulator import SRWUndulator

import scipy.constants as codata
m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)

class FluxPrecisionParameters(object):
    def __init__(self,
                 initial_UR_harmonic = 1,
                 final_UR_harmonic = 21, #Maximum number of harmonics considered. This is critical for speed
                 longitudinal_integration_precision_parameter = 1.5,
                 azimuthal_integration_precision_parameter = 1.5,
                 calculation_type = 2 #calculate flux (1) or flux per unit surface (2)
                 ):
        self._initial_UR_harmonic = initial_UR_harmonic
        self._final_UR_harmonic = final_UR_harmonic
        self._longitudinal_integration_precision_parameter = longitudinal_integration_precision_parameter
        self._azimuthal_integration_precision_parameter = azimuthal_integration_precision_parameter
        self._calculation_type = calculation_type

    def to_SRW_array(self):
        return [int(self._initial_UR_harmonic),
                int(self._final_UR_harmonic),
                float(self._longitudinal_integration_precision_parameter),
                float(self._azimuthal_integration_precision_parameter),
                int(self._calculation_type)]

class SRWUndulatorLightSource(SRWLightSource):

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
                 K_vertical = 0.0,
                 K_horizontal = 0.0,
                 period_length = 0.0,
                 number_of_periods = 1):

        super().__init__(name,
                         electron_energy_in_GeV = electron_energy_in_GeV,
                         electron_energy_spread = electron_energy_spread,
                         ring_current = ring_current,
                         electrons_per_bunch = electrons_per_bunch,
                         electron_beam_size_h=electron_beam_size_h,
                         electron_beam_size_v=electron_beam_size_v,
                         emittance=emittance,
                         coupling_costant=coupling_costant,
                         magnetic_structure=SRWUndulator(K_vertical,
                                                         K_horizontal,
                                                         period_length,
                                                         number_of_periods))

    def get_length(self):
        return self._magnetic_structure._period_length*self._magnetic_structure._number_of_periods

    def get_resonance_wavelength(self):
        return (1 + (self._magnetic_structure._K_vertical**2 + self._magnetic_structure._K_horizontal**2) / 2.0) / 2 / self.get_gamma()**2 * self._magnetic_structure._period_length

    def get_resonance_energy(self):
        return m2ev / self.get_resonance_wavelength()

    def get_photon_source_properties(self, harmonic):
        wavelength = m2ev/(harmonic*self.get_resonance_energy())
        undulator_length = self.get_length()

        # calculate sizes of the photon undulator beam
        # see formulas 25 & 30 in Elleaume (Onaki & Elleaume)
        s_phot = 2.740/(4e0*numpy.pi)*numpy.sqrt(undulator_length*wavelength)
        sp_phot = 0.69*numpy.sqrt(wavelength/undulator_length)

        photon_h = numpy.sqrt(numpy.power(numpy.sqrt(self._electron_beam._moment_xx), 2) + numpy.power(s_phot, 2))
        photon_v = numpy.sqrt(numpy.power(numpy.sqrt(self._electron_beam._moment_yy), 2) + numpy.power(s_phot, 2))
        photon_hp = numpy.sqrt(numpy.power(numpy.sqrt(self._electron_beam._moment_xpxp), 2) + numpy.power(sp_phot, 2))
        photon_vp = numpy.sqrt(numpy.power(numpy.sqrt(self._electron_beam._moment_ypyp), 2) + numpy.power(sp_phot, 2))

        cohH = wavelength/4/numpy.pi / photon_h / photon_hp
        cohV = wavelength/4/numpy.pi / photon_v / photon_vp

        dls = numpy.sqrt(2*undulator_length*wavelength)/4/numpy.pi

        return PhotonSourceProperties(rms_h=photon_h,
                                      rms_v=photon_v,
                                      rms_hp=photon_hp,
                                      rms_vp=photon_vp,
                                      coherence_volume_h=cohH,
                                      coherence_volume_v=cohV,
                                      diffraction_limit=dls)

    def get_undulator_flux(self,
                           source_wavefront_parameters = SourceWavefrontParameters(),
                           flux_precision_parameters = FluxPrecisionParameters()):

        stkF = source_wavefront_parameters.to_SRWLStokes()

        srwl.CalcStokesUR(stkF,
                          self._electron_beam.to_SRWLPartBeam(),
                          self._magnetic_structure.get_SRWMagneticStructure(),
                          flux_precision_parameters.to_SRW_array())

        eArray = numpy.zeros(stkF.mesh.ne)
        intensArray = numpy.zeros(stkF.mesh.ne)

        for i in range(source_wavefront_parameters._photon_energy_points):
            eArray[i] = stkF.mesh.eStart+i*(stkF.mesh.eFin-stkF.mesh.eStart)/numpy.array((stkF.mesh.ne-1)).clip(min=1)
            intensArray[i] = stkF.arS[i]

        return (eArray, intensArray)

