from wofrysrw.storage_ring.srw_light_source import SRWLightSource
from wofrysrw.storage_ring.magnetic_structures.srw_bending_magnet import SRWBendingMagnet

class SRWBendingMagnetLightSource(SRWLightSource):

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
                 magnetic_radius = 0.0,
                 magnetic_field = 0.0,
                 length = 0.0):

        magnetic_radius = abs(magnetic_radius)
        magnetic_field = abs(magnetic_field)

        if magnetic_radius ==0 and magnetic_field ==0: raise ValueError("Magnetic Radius and Magnetic Field cannot be both equal to 0.0")

        magnetic_radius = self.calculate_magnetic_radius(magnetic_field, electron_energy_in_GeV) if magnetic_radius == 0.0 else magnetic_radius
        magnetic_field = self.calculate_magnetic_field(magnetic_radius, electron_energy_in_GeV) if magnetic_field == 0.0 else magnetic_field

        super().__init__(name,
                         electron_energy_in_GeV = electron_energy_in_GeV,
                         electron_energy_spread = electron_energy_spread,
                         ring_current = ring_current,
                         electrons_per_bunch = electrons_per_bunch,
                         electron_beam_size_h=electron_beam_size_h,
                         electron_beam_size_v=electron_beam_size_v,
                         emittance=emittance,
                         coupling_costant=coupling_costant,
                         magnetic_structure=SRWBendingMagnet(magnetic_radius,
                                                             magnetic_field,
                                                             length))

    def calculate_magnetic_field(self, magnetic_radius, electron_energy):
           return 3.334728*electron_energy*1e-3/magnetic_radius

    def calculate_magnetic_radius(self, magnetic_field, electron_energy):
           return 3.334728*electron_energy*1e-3/magnetic_field
