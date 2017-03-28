import numpy

from wofrysrw.storage_ring.srw_light_source import SRWLightSource
from wofrysrw.storage_ring.light_sources.srw_wiggler_light_source import SRWWigglerLightSource
from wofrysrw.storage_ring.srw_light_source import SourceWavefrontParameters, SRWPrecisionParameters

from srxraylib.plot.gol import plot, plot_contour, plot_surface

period_length=0.14

wiggler = SRWWigglerLightSource("XRD1",
                                electron_energy_in_GeV=2.0,
                                electron_energy_spread=0.0007,
                                ring_current=0.4,
                                electron_beam_size_h=0.05545e-3,
                                electron_beam_size_v=2.784e-6,
                                emittance=0.2525e-9,
                                coupling_costant=0.01,
                                K_vertical=21.007,
                                period_length=period_length,
                                number_of_periods=int(4.0/period_length))

print(wiggler.get_electron_beam().get_electron_beam_geometrical_properties().to_info())


wf_parameters = SourceWavefrontParameters(photon_energy_min = 13000,
                                          photon_energy_max = 13000,
                                          photon_energy_points=1,
                                          h_slit_gap = 1.5e-3,
                                          v_slit_gap = 0.2e-3,
                                          h_slit_points=100,
                                          v_slit_points=10,
                                          distance = 2.0,
                                          srw_precision_parameters=SRWPrecisionParameters(sr_method=2,
                                                                                          relative_precision=0.01,
                                                                                          sampling_factor_for_adjusting_nx_ny=0.2))

e, h, v, i = wiggler.get_radiation(wiggler.get_SRW_Wavefront(wf_parameters))


plot_contour(i[int(int(e.size/2))],h*1e3,v*1e3,title="%s SRW; E=%g eV"%("XRD1",e[int(e.size/2)]),xtitle="H [mm]",ytitle="V [mm]",plot_points=0,
             contour_levels=numpy.linspace(0, i.max(), 20), cmap=None, cbar=1,cbar_title="Flux ",show=1)

plot_surface(i[int(e.size/2)],h*1e3,v*1e3,title="%s SRW; E=%g eV"%("XRD1",e[int(e.size/2)]),xtitle="H [mm]",ytitle="V [mm]",show=1)
