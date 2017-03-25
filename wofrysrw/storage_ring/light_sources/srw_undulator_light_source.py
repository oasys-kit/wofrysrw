import numpy

from srwlib import srwl, array, SRWLMagFldU, SRWLMagFldH, SRWLMagFldC, SRWLPartBeam, SRWLStokes

from wofrysrw.storage_ring.srw_light_source import SRWLightSource, PhotonSourceProperties, SourceWavefrontParameters
from wofrysrw.storage_ring.magnetic_structures.srw_undulator import SRWUndulator
from wofrysrw.storage_ring.srw_electron_beam import SRWElectronBeam, SRWElectronBeamGeometricalProperties

import scipy.constants as codata
codata_mee = (codata.m_e * codata.c**2 / codata.e) * 1e-6 # electron mass energy equivalent in MeV
m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)
cte = codata.e/(2*numpy.pi*codata.electron_mass*codata.c)

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

        electron_beam = SRWElectronBeam(energy_in_GeV=electron_energy_in_GeV,
                                        energy_spread=electron_energy_spread,
                                        current=ring_current,
                                        electrons_per_bunch=electrons_per_bunch)

        electron_beam.set_moments_from_electron_beam_geometrical_properties(SRWElectronBeamGeometricalProperties(electron_beam_size_h=electron_beam_size_h,
                                                                                                                 electron_beam_divergence_h=(emittance/electron_beam_size_h),
                                                                                                                 electron_beam_size_v=electron_beam_size_v,
                                                                                                                 electron_beam_divergence_v=(coupling_costant*emittance/electron_beam_size_v)))
        super().__init__(name,
                         electron_beam=electron_beam,
                         magnetic_structure=SRWUndulator(K_vertical,
                                                         K_horizontal,
                                                         period_length,
                                                         number_of_periods))

    def get_length(self):
        return self._magnetic_structure._period_length*self._magnetic_structure._number_of_periods

    def get_gamma(self):
        return self._electron_beam._energy_in_GeV / (codata_mee * 1e-3)

    def get_resonance_wavelength(self):
        return (1 + (self._magnetic_structure._K_vertical**2 + self._magnetic_structure._K_horizontal**2) / 2.0) / 2 / self.get_gamma()**2 * self._magnetic_structure._period_length

    def get_resonance_energy(self):
        return m2ev / self.get_resonance_wavelength()

    # calculate sizes of the photon undulator beam
    # see formulas 25 & 30 in Elleaume (Onaki & Elleaume)
    def get_photon_source_properties(self, harmonic):
        wavelength = m2ev/(harmonic*self.get_resonance_energy())
        undulator_length = self.get_length()

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

    def get_flux(self,
                 source_wavefront_parameters = SourceWavefrontParameters(),
                 max_harmonic_number = 21): #Maximum number of harmonics considered. This is critical for speed


        B0_v = self._magnetic_structure._K_vertical/self._magnetic_structure._period_length/cte
        B0_h = self._magnetic_structure._K_horizontal/self._magnetic_structure._period_length/cte

        #***********Undulator
        harmB_v = SRWLMagFldH() #magnetic field harmonic
        harmB_v.n = 1 #harmonic number
        harmB_v.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
        harmB_v.B = B0_v #magnetic field amplitude [T]

        harmB_h = SRWLMagFldH() #magnetic field harmonic
        harmB_h.n = 1 #harmonic number
        harmB_h.h_or_v = 'h' #magnetic field plane: horzontal ('h') or vertical ('v')
        harmB_h.B = B0_h #magnetic field amplitude [T]

        und = SRWLMagFldU([harmB_v, harmB_h])
        und.per = self._magnetic_structure._period_length #period length [m]
        und.nPer = self._magnetic_structure._number_of_periods #number of periods (will be rounded to integer)

        #***********Electron Beam
        eBeam = SRWLPartBeam()
        eBeam.Iavg = self._electron_beam._current #average current [A]
        eBeam.partStatMom1.x = 0. #initial transverse positions [m]
        eBeam.partStatMom1.y = 0.
        eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
        eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
        eBeam.partStatMom1.yp = 0
        eBeam.partStatMom1.gamma = self._electron_beam._energy_in_GeV*1e3/codata_mee #relative energy

        #2nd order stat. moments:
        eBeam.arStatMom2[0] = self._electron_beam._moment_xx #<(x-<x>)^2>
        eBeam.arStatMom2[1] = self._electron_beam._moment_xxp #<(x-<x>)(x'-<x'>)>
        eBeam.arStatMom2[2] = self._electron_beam._moment_xpxp #<(x'-<x'>)^2>
        eBeam.arStatMom2[3] = self._electron_beam._moment_yy #<(y-<y>)^2>
        eBeam.arStatMom2[4] = self._electron_beam._moment_yyp #<(y-<y>)(y'-<y'>)>
        eBeam.arStatMom2[5] = self._electron_beam._moment_ypyp #<(y'-<y'>)^2>
        eBeam.arStatMom2[10] = self._electron_beam._energy_spread**2 #<(E-<E>)^2>/<E>^2

        #***********Precision Parameters
        arPrecF = [0]*5 #for spectral flux vs photon energy
        arPrecF[0] = 1 #initial UR harmonic to take into account
        arPrecF[1] = max_harmonic_number #final UR harmonic to take into account
        arPrecF[2] = 1.5 #longitudinal integration precision parameter
        arPrecF[3] = 1.5 #azimuthal integration precision parameter
        arPrecF[4] = 1 #calculate flux (1) or flux per unit surface (2)

        #***********UR Stokes Parameters (mesh) for Spectral Flux
        stkF = SRWLStokes() #for spectral flux vs photon energy
        stkF.allocate(source_wavefront_parameters._photon_energy_points, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
        stkF.mesh.zStart = source_wavefront_parameters._distance #longitudinal position [m] at which UR has to be calculated
        stkF.mesh.eStart = source_wavefront_parameters._photon_energy_min #initial photon energy [eV]
        stkF.mesh.eFin =   source_wavefront_parameters._photon_energy_max #final photon energy [eV]
        stkF.mesh.xStart = -source_wavefront_parameters._h_slit_gap/2 #initial horizontal position [m]
        stkF.mesh.xFin =    source_wavefront_parameters._h_slit_gap/2 #final horizontal position [m]
        stkF.mesh.yStart = -source_wavefront_parameters._v_slit_gap/2 #initial vertical position [m]
        stkF.mesh.yFin =    source_wavefront_parameters._v_slit_gap/2 #final vertical position [m]

        #**********************Calculation (SRWLIB function calls)
        print('Performing Spectral Flux (Stokes parameters) calculation ... ') # , end='')

        srwl.CalcStokesUR(stkF, eBeam, und, arPrecF)

        eArray = numpy.zeros(source_wavefront_parameters._photon_energy_points)
        intensArray = numpy.zeros(source_wavefront_parameters._photon_energy_points)
        for i in range(stkF.mesh.ne):
            ener = stkF.mesh.eStart+i*(stkF.mesh.eFin-stkF.mesh.eStart)/numpy.array((stkF.mesh.ne-1)).clip(min=1)
            eArray[i] = ener
            intensArray[i] = stkF.arS[i]

        return (eArray, intensArray)

    def get_power_density(self,
                 source_wavefront_parameters = SourceWavefrontParameters(),
                 max_harmonic_number = 21): #Maximum number of harmonics considered. This is critical for speed
        global scanCounter

        B0_v = self._magnetic_structure._K_vertical/self._magnetic_structure._period_length/cte
        B0_h = self._magnetic_structure._K_horizontal/self._magnetic_structure._period_length/cte

        #***********Undulator
        harmB_v = SRWLMagFldH() #magnetic field harmonic
        harmB_v.n = 1 #harmonic number
        harmB_v.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
        harmB_v.B = B0_v #magnetic field amplitude [T]

        harmB_h = SRWLMagFldH() #magnetic field harmonic
        harmB_h.n = 1 #harmonic number
        harmB_h.h_or_v = 'h' #magnetic field plane: horzontal ('h') or vertical ('v')
        harmB_h.B = B0_h #magnetic field amplitude [T]

        und = SRWLMagFldU([harmB_v, harmB_h])
        und.per = self._magnetic_structure._period_length #period length [m]
        und.nPer = self._magnetic_structure._number_of_periods #number of periods (will be rounded to integer)

        magFldCnt = SRWLMagFldC([und], array('d', [0]), array('d', [0]), array('d', [0]))

        #***********Electron Beam
        eBeam = SRWLPartBeam()
        eBeam.Iavg = self._electron_beam._current #average current [A]
        eBeam.partStatMom1.x = 0. #initial transverse positions [m]
        eBeam.partStatMom1.y = 0.
        eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
        eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
        eBeam.partStatMom1.yp = 0
        eBeam.partStatMom1.gamma = self._electron_beam._energy_in_GeV*1e3/codata_mee #relative energy

        #2nd order stat. moments:
        eBeam.arStatMom2[0] = self._electron_beam._moment_xx #<(x-<x>)^2>
        eBeam.arStatMom2[1] = self._electron_beam._moment_xxp #<(x-<x>)(x'-<x'>)>
        eBeam.arStatMom2[2] = self._electron_beam._moment_xpxp #<(x'-<x'>)^2>
        eBeam.arStatMom2[3] = self._electron_beam._moment_yy #<(y-<y>)^2>
        eBeam.arStatMom2[4] = self._electron_beam._moment_yyp #<(y-<y>)(y'-<y'>)>
        eBeam.arStatMom2[5] = self._electron_beam._moment_ypyp #<(y'-<y'>)^2>
        eBeam.arStatMom2[10] = self._electron_beam._energy_spread**2 #<(E-<E>)^2>/<E>^2

        #***********Precision Parameters
        arPrecP = [0]*5 #for power density
        arPrecP[0] = 1.5 #precision factor
        arPrecP[1] = 1 #power density computation method (1- "near field", 2- "far field")
        arPrecP[2] = 0 #initial longitudinal position (effective if arPrecP[2] < arPrecP[3])
        arPrecP[3] = 0 #final longitudinal position (effective if arPrecP[2] < arPrecP[3])
        arPrecP[4] = 20000 #number of points for (intermediate) trajectory calculation

        #***********UR Stokes Parameters (mesh) for power densiyu
        stkP = SRWLStokes() #for power density
        stkP.allocate(source_wavefront_parameters._photon_energy_points, source_wavefront_parameters._h_slit_points,source_wavefront_parameters._v_slit_points) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
        stkP.mesh.eStart = source_wavefront_parameters._photon_energy_min #initial photon energy [eV]
        stkP.mesh.eFin =   source_wavefront_parameters._photon_energy_max #final photon energy [eV]
        stkP.mesh.zStart = source_wavefront_parameters._distance #longitudinal position [m] at which power density has to be calculated
        stkP.mesh.xStart = -source_wavefront_parameters._h_slit_gap/2 #initial horizontal position [m]
        stkP.mesh.xFin =    source_wavefront_parameters._h_slit_gap/2 #final horizontal position [m]
        stkP.mesh.yStart = -source_wavefront_parameters._v_slit_gap/2 #initial vertical position [m]
        stkP.mesh.yFin =    source_wavefront_parameters._v_slit_gap/2 #final vertical position [m]

        #**********************Calculation (SRWLIB function calls)
        srwl.CalcPowDenSR(stkP, eBeam, 0, magFldCnt, arPrecP)

        hArray = numpy.zeros(stkP.mesh.nx)
        vArray = numpy.zeros(stkP.mesh.ny)
        totPower = numpy.array(0.0)

        powerArray = numpy.zeros((stkP.mesh.nx,stkP.mesh.ny))

        # fill arrays
        ij = -1
        for j in range(stkP.mesh.ny):
            for i in range(stkP.mesh.nx):
                ij += 1
                xx = stkP.mesh.xStart + i*(stkP.mesh.xFin-stkP.mesh.xStart)/(stkP.mesh.nx-1)
                yy = stkP.mesh.yStart + j*(stkP.mesh.yFin-stkP.mesh.yStart)/(stkP.mesh.ny-1)
                totPower += stkP.arS[ij]
                powerArray[i,j] = stkP.arS[ij]
                hArray[i] = xx*1e3 # mm
                vArray[j] = yy*1e3 # mm

        return (hArray, vArray, powerArray)
