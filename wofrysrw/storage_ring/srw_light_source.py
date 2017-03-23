import numpy
from array import array as _array
import scipy.constants as codata

from srwlib import *

from syned.storage_ring.light_source import LightSource

from wofrysrw.storage_ring.srw_magnetic_structure import SRWMagneticStructure
from wofrysrw.storage_ring.srw_electron_beam import SRWElectronBeam
from wofrysrw.propagator.wavefront2D.srw_wavefront import SRWWavefront

class SourceWavefrontParameters(object):
    def __init__(self, 
                 photon_energy_min = 100,
                 photon_energy_max = 10100,
                 photon_energy_points = 51,
                 h_slit_gap = 1e-3, 
                 h_slit_points = 51, 
                 v_slit_gap = 1e-3, 
                 v_slit_points = 51, 
                 distance = 10):
        self._photon_energy_min    = photon_energy_min    
        self._photon_energy_max    = photon_energy_max   
        self._photon_energy_points = photon_energy_points
        self._h_slit_gap           = h_slit_gap          
        self._h_slit_points        = h_slit_points       
        self._v_slit_gap           = v_slit_gap          
        self._v_slit_points        = v_slit_points        
        self._distance             = distance            

class PhotonSourceProperties(object):

    def __init__(self,
                 rms_h = 0.0,
                 rms_v = 0.0,
                 rms_hp = 0.0,
                 rms_vp = 0.0,
                 coherence_volume_h = 0.0,
                 coherence_volume_v = 0.0,
                 diffraction_limit = 0.0):

        self._rms_h = rms_h,
        self._rms_v = rms_v,
        self._rms_hp = rms_hp,
        self._rms_vp = rms_vp,
        self._coherence_volume_h = coherence_volume_h
        self._coherence_volume_v = coherence_volume_v
        self._diffraction_limit = diffraction_limit

class SRWLightSource(LightSource):
    def __init__(self, name="Undefined", electron_beam=SRWElectronBeam(), magnetic_structure=SRWMagneticStructure()):
        LightSource.__init__(self, name, electron_beam, magnetic_structure)

    def get_photon_source_properties(self):
        return NotImplementedError("must be implemented in subclasses")

    def get_SRW_Wavefront(self, source_wavefront_parameters = SourceWavefrontParameters()):

        mesh = SRWLRadMesh(source_wavefront_parameters._photon_energy_min,
                           source_wavefront_parameters._photon_energy_max,
                           source_wavefront_parameters._photon_energy_points,
                           -source_wavefront_parameters._h_slit_gap/2,source_wavefront_parameters._h_slit_gap/2,source_wavefront_parameters._h_slit_points,
                           -source_wavefront_parameters._v_slit_gap/2,source_wavefront_parameters._v_slit_gap/2,source_wavefront_parameters._v_slit_points, 
                           source_wavefront_parameters._distance)


        paramSE = [1, 0.01, 0, 0, 50000, 1, 0]

        wfr = SRWLWfr()
        wfr.mesh = mesh
        wfr.partBeam = self._electron_beam.to_SRWLPartBeam()
        wfr.allocate(mesh.ne, mesh.nx, mesh.ny)

        srwl.CalcElecFieldSR(wfr, 0, self._magnetic_structure.get_SRWLMagFldC(), paramSE)

        return SRWWavefront.decorateSRWWF(wfr)

    def get_intensity(self, srw_wavefront):

        mesh0 = srw_wavefront.mesh

        INTENSITY_TYPE_MULTI_ELECTRON=1

        hArray=numpy.linspace(srw_wavefront.mesh.xStart,srw_wavefront.mesh.xFin, srw_wavefront.mesh.nx)
        vArray=numpy.linspace(srw_wavefront.mesh.yStart,srw_wavefront.mesh.yFin, srw_wavefront.mesh.ny)
        eArray=numpy.linspace(srw_wavefront.mesh.eStart,srw_wavefront.mesh.eFin, srw_wavefront.mesh.ne)

        intensArray = numpy.zeros((eArray.size,hArray.size,vArray.size,))
        for ie in range(eArray.size):
            arI0 = _array('f', [0]*mesh0.nx*mesh0.ny) #"flat" array to take 2D intensity data
            srwl.CalcIntFromElecField(arI0, srw_wavefront, 6, INTENSITY_TYPE_MULTI_ELECTRON, 3, eArray[ie], 0, 0) # 6 is for total polarizarion; 0=H, 1=V

            data = numpy.ndarray(buffer=arI0, shape=(mesh0.ny, mesh0.nx),dtype=arI0.typecode)

            for ix in range(hArray.size):
                for iy in range(vArray.size):
                    intensArray[ie,ix,iy,] = data[iy,ix]
        
        return (eArray, hArray, vArray, intensArray)

    @classmethod
    def get_spectrum_from_intensity(cls, hArray, vArray, intensArray):
        return (intensArray.sum(axis=2)).sum(axis=1)*(vArray[1]-vArray[0])*(hArray[1]-hArray[0])

    @classmethod
    def get_power_density_from_intensity(cls, eArray, intensArray):
        return intensArray.sum(axis=0)*(eArray[1]-eArray[0])*codata.e*1e3
