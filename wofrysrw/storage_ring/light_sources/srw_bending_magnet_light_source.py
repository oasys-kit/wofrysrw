from wofrysrw.storage_ring.srw_light_source import SRWLightSource
from wofrysrw.storage_ring.magnetic_structures.srw_bending_magnet import SRWBendingMagnet

from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet

class SRWBendingMagnetLightSource(SRWLightSource):

    def __init__(self,
                 name="Undefined",
                 electron_energy_in_GeV = 1.0,
                 electron_energy_spread = 0.0,
                 ring_current = 0.1,
                 number_of_bunches = 400,
                 electron_beam_size_h=0.0,
                 electron_beam_size_v=0.0,
                 electron_beam_divergence_h=0.0,
                 electron_beam_divergence_v=0.0,
                 magnetic_radius = 0.0,
                 magnetic_field = 0.0,
                 length = 0.0):

        magnetic_radius = abs(magnetic_radius)
        magnetic_field = abs(magnetic_field)

        if magnetic_radius ==0 and magnetic_field ==0: raise ValueError("Magnetic Radius and Magnetic Field cannot be both equal to 0.0")

        magnetic_radius = BendingMagnet.calculate_magnetic_radius(magnetic_field, electron_energy_in_GeV) if magnetic_radius == 0.0 else magnetic_radius
        magnetic_field = BendingMagnet.calculate_magnetic_field(magnetic_radius, electron_energy_in_GeV) if magnetic_field == 0.0 else magnetic_field

        super().__init__(name,
                         electron_energy_in_GeV = electron_energy_in_GeV,
                         electron_energy_spread = electron_energy_spread,
                         ring_current = ring_current,
                         number_of_bunches = number_of_bunches,
                         electron_beam_size_h=electron_beam_size_h,
                         electron_beam_size_v=electron_beam_size_v,
                         electron_beam_divergence_h=electron_beam_divergence_h,
                         electron_beam_divergence_v=electron_beam_divergence_v,
                         magnetic_structure=SRWBendingMagnet(magnetic_radius,
                                                             magnetic_field,
                                                             length))
