
from wofrysrw.storage_ring.srw_light_source import SRWLightSource
from wofrysrw.storage_ring.light_sources.srw_undulator_light_source import SRWUndulatorLightSource
from wofrysrw.storage_ring.srw_light_source import SourceWavefrontParameters

from srxraylib.plot.gol import plot, plot_contour

period_length=0.02
total_length = 1.5

print(1)

undulator = SRWUndulatorLightSource("MicroXRD",
                                    electron_energy_in_GeV=2.0,
                                    electron_energy_spread=0.0007,
                                    ring_current=0.4,
                                    electron_beam_size_h=0.005545e-3,
                                    electron_beam_size_v=2.784e-6,
                                    emittance=0.2525e-9,
                                    coupling_costant=0.01,
                                    K_vertical=1.5,
                                    period_length=period_length,
                                    number_of_periods=int(total_length/period_length))

print(2)

wf_parameters = SourceWavefrontParameters(photon_energy_min = 1,
                                          photon_energy_max = 12001,
                                          photon_energy_points=201,
                                          h_slit_gap = 1e-3,
                                          v_slit_gap = 1e-3,
                                          distance = 10.0)

print(3)

wavefront = undulator.get_SRW_Wavefront(wf_parameters)

print(4)

e, h, v, i = undulator.get_intensity(wavefront)

print(5)

plot(e, SRWLightSource.get_spectrum_from_intensity(h, v, i), show=1, title="Spectrum for MicroXRD")
plot_contour(SRWLightSource.get_power_density_from_intensity(e, i), h, v, title="Power density", show=1)
