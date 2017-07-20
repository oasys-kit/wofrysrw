from syned.beamline.shape import BoundaryShape

from wofrysrw.beamline.optical_elements.absorbers.srw_slit import SRWSlit

class SRWAperture(SRWSlit):
    def __init__(self, name="Undefined", boundary_shape=BoundaryShape()):
        SRWSlit.__init__(self, name=name, boundary_shape=boundary_shape)

    def get_srw_ap_or_ob(self):
       return 'a'