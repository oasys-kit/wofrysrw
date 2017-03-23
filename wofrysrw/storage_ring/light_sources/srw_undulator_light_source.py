import numpy

from wofrysrw.storage_ring.srw_light_source import SRWLightSource, PhotonSourceProperties
from wofrysrw.storage_ring.magnetic_structures.srw_undulator import SRWUndulator
from wofrysrw.storage_ring.srw_electron_beam import SRWElectronBeam

import scipy.constants as codata
codata_mee = (codata.m_e * codata.c**2 / codata.e) * 1e-6 # electron mass energy equivalent in MeV
m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)

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
                         electron_beam=SRWElectronBeam(electron_energy_in_GeV,
                                                       electron_energy_spread,
                                                       ring_current,
                                                       electrons_per_bunch,
                                                       electron_beam_size_h**2,
                                                       0.0,
                                                       (emittance/electron_beam_size_h)**2,
                                                       electron_beam_size_v**2,
                                                       0.0,
                                                       (coupling_costant*emittance/electron_beam_size_v)**2),
                         magnetic_structure=SRWUndulator(K_vertical,
                                                         K_horizontal,
                                                         period_length,
                                                         number_of_periods))

    def get_length(self):
        return self._magnetic_structure._period_length*self._magnetic_structure.number_of_periods

    def get_gamma(self):
        self._electron_beam._energy_in_GeV / (codata_mee * 1e-3)

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

    def get_flux(self, photon_energy_min = 100,
                 photon_energy_max = 100100,
                 photon_energy_points = 1000, max_harmonic_number = 21):
        pass
    '''
    global scanCounter

    t0 = time.time()
    print("Inside calc1d_srw")
    #Maximum number of harmonics considered. This is critical for speed.
    #TODO: set it automatically to a reasonable value (see how is done by Urgent).
    Nmax = srw_max_harmonic_number # 21,61
    #derived
    #TODO calculate the numerical factor using codata
    #B0 = bl['Kv']/0.934/(bl['PeriodID']*1e2)

    cte = codata.e/(2*numpy.pi*codata.electron_mass*codata.c)
    B0 = bl['Kv']/bl['PeriodID']/cte


    print('Running SRW (SRWLIB Python)')

    #***********Undulator
    harmB = srwlib.SRWLMagFldH() #magnetic field harmonic
    harmB.n = 1 #harmonic number
    harmB.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
    harmB.B = B0 #magnetic field amplitude [T]

    und = srwlib.SRWLMagFldU([harmB])
    und.per = bl['PeriodID'] #period length [m]
    und.nPer = bl['NPeriods'] #number of periods (will be rounded to integer)

    #Container of all magnetic field elements
    magFldCnt = srwlib.SRWLMagFldC([und], srwlib.array('d', [0]), srwlib.array('d', [0]), srwlib.array('d', [0]))

    #***********Electron Beam
    eBeam = srwlib.SRWLPartBeam()
    eBeam.Iavg = bl['ElectronCurrent'] #average current [A]
    eBeam.partStatMom1.x = 0. #initial transverse positions [m]
    eBeam.partStatMom1.y = 0.
    eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
    eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
    eBeam.partStatMom1.yp = 0
    eBeam.partStatMom1.gamma = bl['ElectronEnergy']*1e3/codata_mee #relative energy

    if zero_emittance:
        sigX     = 1e-25
        sigXp    = 1e-25
        sigY     = 1e-25
        sigYp    = 1e-25
        sigEperE = 1e-25
    else:
        sigX =  bl['ElectronBeamSizeH'] #horizontal RMS size of e-beam [m]
        sigXp = bl['ElectronBeamDivergenceH'] #horizontal RMS angular divergence [rad]
        sigY =  bl['ElectronBeamSizeV'] #vertical RMS size of e-beam [m]
        sigYp = bl['ElectronBeamDivergenceV'] #vertical RMS angular divergence [rad]
        sigEperE = bl['ElectronEnergySpread']

    print("calc1dSrw: starting calculation using ElectronEnergySpead=%e \n"%((sigEperE)))

    #2nd order stat. moments:
    eBeam.arStatMom2[0] = sigX*sigX #<(x-<x>)^2>
    eBeam.arStatMom2[1] = 0 #<(x-<x>)(x'-<x'>)>
    eBeam.arStatMom2[2] = sigXp*sigXp #<(x'-<x'>)^2>
    eBeam.arStatMom2[3] = sigY*sigY #<(y-<y>)^2>
    eBeam.arStatMom2[4] = 0 #<(y-<y>)(y'-<y'>)>
    eBeam.arStatMom2[5] = sigYp*sigYp #<(y'-<y'>)^2>
    eBeam.arStatMom2[10] = sigEperE*sigEperE #<(E-<E>)^2>/<E>^2

    #***********Precision Parameters
    arPrecF = [0]*5 #for spectral flux vs photon energy
    arPrecF[0] = 1 #initial UR harmonic to take into account
    arPrecF[1] = Nmax #final UR harmonic to take into account
    arPrecF[2] = 1.5 #longitudinal integration precision parameter
    arPrecF[3] = 1.5 #azimuthal integration precision parameter
    arPrecF[4] = 1 #calculate flux (1) or flux per unit surface (2)

    #***********UR Stokes Parameters (mesh) for Spectral Flux
    stkF = srwlib.SRWLStokes() #for spectral flux vs photon energy
    #srio stkF.allocate(10000, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
    stkF.allocate(photonEnergyPoints, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
    stkF.mesh.zStart = bl['distance'] #longitudinal position [m] at which UR has to be calculated
    stkF.mesh.eStart = photonEnergyMin #initial photon energy [eV]
    stkF.mesh.eFin =   photonEnergyMax #final photon energy [eV]
    stkF.mesh.xStart = -bl['gapH']/2 #initial horizontal position [m]
    stkF.mesh.xFin =    bl['gapH']/2 #final horizontal position [m]
    stkF.mesh.yStart = -bl['gapV']/2 #initial vertical position [m]
    stkF.mesh.yFin =    bl['gapV']/2 #final vertical position [m]

    #**********************Calculation (SRWLIB function calls)
    print('Performing Spectral Flux (Stokes parameters) calculation ... ') # , end='')

    srwlib.srwl.CalcStokesUR(stkF, eBeam, und, arPrecF)

    print('Done calc1dSrw calculation in sec '+str(time.time()-t0))
    #**********************Saving results

    if fileName is not None:
        if fileAppend:
            f = open(fileName,"a")
        else:
            scanCounter = 0
            f = open(fileName,"w")
            f.write("#F "+fileName+"\n")

        f.write("\n")
        scanCounter +=1
        f.write("#S %d Undulator spectrum calculation using SRW\n"%(scanCounter))

        for i,j in bl.items(): # write bl values
            f.write ("#UD %s = %s\n" % (i,j) )
        f.write("#UD photonEnergyMin =  %f\n"%(photonEnergyMin))
        f.write("#UD photonEnergyMax =  %f\n"%(photonEnergyMax))
        f.write("#UD photonEnergyPoints =  %d\n"%(photonEnergyPoints))
        f.write("#UD B0 =  %f\n"%(B0))

        #
        # write flux to file
        #
        header="#N 4 \n#L PhotonEnergy[eV]  PhotonWavelength[A]  Flux[phot/sec/0.1%bw]  Spectral Power[W/eV]\n"
        f.write(header)

    eArray = numpy.zeros(photonEnergyPoints)
    intensArray = numpy.zeros(photonEnergyPoints)
    for i in range(stkF.mesh.ne):
        ener = stkF.mesh.eStart+i*(stkF.mesh.eFin-stkF.mesh.eStart)/numpy.array((stkF.mesh.ne-1)).clip(min=1)
        if fileName is not None: f.write(' ' + repr(ener) + '   ' + repr(m2ev/ener*1e10) + '    ' +
                repr(stkF.arS[i]) + '    ' +
                repr(stkF.arS[i]*codata.e*1e3) + '\n')
        eArray[i] = ener
        intensArray[i] = stkF.arS[i]

    if fileName is not None:
        f.close()

        if fileAppend:
            print("Data appended to file: "+fileName)
        else:
            print("File written to disk: "+fileName)


    return (eArray,intensArray)
    '''

    def get_power_density(self):
        pass
    '''
    global scanCounter
    print("Inside calc2d_srw")
    #Maximum number of harmonics considered. This is critical for speed.
    #TODO: set it automatically to a reasonable value (see how is done by Urgent).
    Nmax = srw_max_harmonic_number # 21,61
    #derived
    #TODO calculate the numerical factor using codata
    # B0 = bl['Kv']/0.934/(bl['PeriodID']*1e2)
    cte = codata.e/(2*numpy.pi*codata.electron_mass*codata.c)
    B0 = bl['Kv']/bl['PeriodID']/cte

    print('Running SRW (SRWLIB Python)')

    #***********Undulator
    harmB = srwlib.SRWLMagFldH() #magnetic field harmonic
    harmB.n = 1 #harmonic number
    harmB.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
    harmB.B = B0 #magnetic field amplitude [T]

    und = srwlib.SRWLMagFldU([harmB])
    und.per = bl['PeriodID'] #period length [m]
    und.nPer = bl['NPeriods'] #number of periods (will be rounded to integer)

    #Container of all magnetic field elements
    magFldCnt = srwlib.SRWLMagFldC([und], array.array('d', [0]), array.array('d', [0]), array.array('d', [0]))

    #***********Electron Beam
    eBeam = srwlib.SRWLPartBeam()
    eBeam.Iavg = bl['ElectronCurrent'] #average current [A]
    eBeam.partStatMom1.x = 0. #initial transverse positions [m]
    eBeam.partStatMom1.y = 0.
    eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
    eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
    eBeam.partStatMom1.yp = 0
    eBeam.partStatMom1.gamma = bl['ElectronEnergy']*1e3/codata_mee #relative energy

    if zero_emittance:
        sigEperE = 1e-25
        sigX     = 1e-25
        sigXp    = 1e-25
        sigY     = 1e-25
        sigYp    = 1e-25
    else:
        sigEperE = bl['ElectronEnergySpread'] #relative RMS energy spread
        sigX =  bl['ElectronBeamSizeH'] #horizontal RMS size of e-beam [m]
        sigXp = bl['ElectronBeamDivergenceH'] #horizontal RMS angular divergence [rad]
        sigY =  bl['ElectronBeamSizeV'] #vertical RMS size of e-beam [m]
        sigYp = bl['ElectronBeamDivergenceV'] #vertical RMS angular divergence [rad]

    #2nd order stat. moments:
    eBeam.arStatMom2[0] = sigX*sigX #<(x-<x>)^2>
    eBeam.arStatMom2[1] = 0 #<(x-<x>)(x'-<x'>)>
    eBeam.arStatMom2[2] = sigXp*sigXp #<(x'-<x'>)^2>
    eBeam.arStatMom2[3] = sigY*sigY #<(y-<y>)^2>
    eBeam.arStatMom2[4] = 0 #<(y-<y>)(y'-<y'>)>
    eBeam.arStatMom2[5] = sigYp*sigYp #<(y'-<y'>)^2>
    eBeam.arStatMom2[10] = sigEperE*sigEperE #<(E-<E>)^2>/<E>^2

    #***********Precision Parameters
    arPrecP = [0]*5 #for power density
    arPrecP[0] = 1.5 #precision factor
    arPrecP[1] = 1 #power density computation method (1- "near field", 2- "far field")
    arPrecP[2] = 0 #initial longitudinal position (effective if arPrecP[2] < arPrecP[3])
    arPrecP[3] = 0 #final longitudinal position (effective if arPrecP[2] < arPrecP[3])
    arPrecP[4] = 20000 #number of points for (intermediate) trajectory calculation

    #***********UR Stokes Parameters (mesh) for power densiyu
    stkP = srwlib.SRWLStokes() #for power density
    stkP.allocate(1, hSlitPoints, vSlitPoints) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
    stkP.mesh.zStart = bl['distance'] #longitudinal position [m] at which power density has to be calculated
    stkP.mesh.xStart = -bl['gapH']/2 #initial horizontal position [m]
    stkP.mesh.xFin =    bl['gapH']/2 #final horizontal position [m]
    stkP.mesh.yStart = -bl['gapV']/2 #initial vertical position [m]
    stkP.mesh.yFin =    bl['gapV']/2 #final vertical position [m]

    #**********************Calculation (SRWLIB function calls)
    print('Performing Power Density calculation (from field) ... ')
    t0 = time.time()
    srwlib.srwl.CalcPowDenSR(stkP, eBeam, 0, magFldCnt, arPrecP)
    print('Done Performing Power Density calculation (from field).')

    #**********************Saving results

    if fileName is not None:
        if fileAppend:
            f = open(fileName,"a")
        else:
            scanCounter = 0
            f = open(fileName,"w")
            f.write("#F "+fileName+"\n")

        #
        # write power density to file as mesh scan
        #
        scanCounter +=1
        f.write("\n#S %d Undulator power density calculation using SRW\n"%(scanCounter))
        for i,j in bl.items(): # write bl values
            f.write ("#UD %s = %s\n" % (i,j) )
        f.write('\n#U B0 = ' + repr(B0 ) + '\n' )
        f.write('\n#U hSlitPoints = ' + repr(hSlitPoints) + '\n' )
        f.write('\n#U vSlitPoints = ' + repr(vSlitPoints) + '\n' )
        f.write("#N 3 \n#L H[mm]  V[mm]  PowerDensity[W/mm^2] \n" )

    hArray = numpy.zeros(stkP.mesh.nx)
    vArray = numpy.zeros(stkP.mesh.ny)
    totPower = numpy.array(0.0)

    hProfile = numpy.zeros(stkP.mesh.nx)
    vProfile = numpy.zeros(stkP.mesh.ny)
    powerArray = numpy.zeros((stkP.mesh.nx,stkP.mesh.ny))

    # fill arrays
    ij = -1
    for j in range(stkP.mesh.ny):
        for i in range(stkP.mesh.nx):
            ij += 1
            xx = stkP.mesh.xStart + i*(stkP.mesh.xFin-stkP.mesh.xStart)/(stkP.mesh.nx-1)
            yy = stkP.mesh.yStart + j*(stkP.mesh.yFin-stkP.mesh.yStart)/(stkP.mesh.ny-1)
            #ij = i*stkP.mesh.nx + j
            totPower += stkP.arS[ij]
            powerArray[i,j] = stkP.arS[ij]
            hArray[i] = xx*1e3 # mm
            vArray[j] = yy*1e3 # mm

    # dump
    if fileName is not None:
        for i in range(stkP.mesh.nx):
            for j in range(stkP.mesh.ny):
                f.write(repr(hArray[i]) + ' ' + repr(vArray[j]) + ' ' + repr(powerArray[i,j]) + '\n')


    totPower = totPower * \
               (stkP.mesh.xFin-stkP.mesh.xStart)/(stkP.mesh.nx-1)*1e3 * \
               (stkP.mesh.yFin-stkP.mesh.yStart)/(stkP.mesh.ny-1)*1e3
    hStep = (stkP.mesh.xFin-stkP.mesh.xStart)/(stkP.mesh.nx-1)

    # dump profiles
    if fileName is not None:


        scanCounter +=1
        f.write("\n#S %d Undulator power density calculation using SRW: H profile\n"%(scanCounter))
        for i,j in bl.items(): # write bl values
            f.write ("#UD %s = %s\n" % (i,j) )
        f.write( "#UD Total power [W]: "+repr(totPower)+"\n")
        f.write( "#UD FWHM [mm] : "+repr(calc_fwhm(hProfile,hStep)[0]*1e3)+"\n")
        f.write( "#N 2 \n")
        f.write( "#L H[mm]  PowerDensityCentralProfile[W/mm2] \n" )
        for i in range(stkP.mesh.nx):
            #xx = stkP.mesh.xStart + i*hStep
            #f.write(repr(xx*1e3) + ' ' + repr(hProfile[i]) + '\n')
            f.write(repr(hArray[i]) + ' ' + \
                    repr(powerArray[i,int(len(vArray)/2)]) + '\n')

        scanCounter +=1
        vStep = (stkP.mesh.yFin-stkP.mesh.yStart)/(stkP.mesh.ny-1)
        f.write("\n#S %d Undulator power density calculation using SRW: V profile\n"%(scanCounter))
        for i,j in bl.items(): # write bl values
            f.write ("#UD %s = %s\n" % (i,j) )
        f.write( "#UD Total power [W]: "+repr(totPower)+"\n")
        f.write( "#UD FWHM [mm] : "+repr(calc_fwhm(vProfile,vStep)[0]*1e3)+"\n")
        f.write( "#N 2 \n")
        f.write( "#L V[mm]  PowerDensityCentralProfile[W/mm2] \n" )
        for j in range(stkP.mesh.ny):
            f.write(repr(vArray[j]) + ' ' +  \
                    repr(powerArray[int(len(hArray)/2),j]) + '\n')

        f.close()

        if fileAppend:
            print("Data appended to file: "+fileName)
        else:
            print("File written to disk: "+fileName)

    print( "Total power SRW [W]: "+repr(totPower))

    return (hArray, vArray, powerArray)
    '''