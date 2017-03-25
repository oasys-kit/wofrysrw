import numpy

from wofrysrw.storage_ring.srw_light_source import SRWLightSource
from wofrysrw.storage_ring.light_sources.srw_undulator_light_source import SRWUndulatorLightSource
from wofrysrw.storage_ring.srw_light_source import SourceWavefrontParameters

from srxraylib.plot.gol import plot, plot_contour, plot_surface

period_length=0.02
total_length = 1.5

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

resonance_energy = undulator.get_resonance_energy()

wf_parameters = SourceWavefrontParameters(photon_energy_min = 1,
                                          photon_energy_max = 12001,
                                          photon_energy_points=101,
                                          h_slit_gap = 5e-3,
                                          v_slit_gap = 5e-3,
                                          distance = 10.0)

wavefront = undulator.get_SRW_Wavefront(wf_parameters)

e, h, v, i = undulator.get_radiation(wavefront)

properties = undulator.get_photon_source_properties(harmonic=1)

print(properties.to_info())

plot(e, SRWLightSource.get_spectrum_from_radiation(h, v, i), show=1, title="Spectrum for MicroXRD")
plot_contour(SRWLightSource.get_power_density_from_radiation(e, i), h, v, title="Power density", show=1)

e, i = undulator.get_flux(source_wavefront_parameters=wf_parameters, max_harmonic_number=7)

plot(e, i, show=1, title="Flux for MicroXRD")

harmonic = 5

wf_parameters._photon_energy_min=harmonic*resonance_energy
wf_parameters._photon_energy_max=harmonic*resonance_energy
wf_parameters._photon_energy_points = 1
wf_parameters._h_slit_gap = 5e-2
wf_parameters._v_slit_gap = 5e-2


h, v, p = undulator.get_power_density(source_wavefront_parameters=wf_parameters, max_harmonic_number=harmonic)


contour_levels = numpy.linspace(0,numpy.max([p.max()]),100)

plot_contour(p,h,v,title="%s SRW, total power = %g"%("Micro XRD", total_power),xtitle="H [mm]",ytitle="V [mm]",plot_points=0,
             contour_levels=contour_levels,cmap=None,cbar=1,cbar_title="Power density [$W/mm^2$]",show=1)

plot_surface(p,h,v,  title="%s SRW, total power = %g"%("Micro XRD", total_power),xtitle="H [mm]",ytitle="V [mm]",show=1)
