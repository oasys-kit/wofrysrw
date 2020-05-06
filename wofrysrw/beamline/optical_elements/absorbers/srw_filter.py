import numpy, xraylib
from scipy.interpolate import RectBivariateSpline
from oasys.util.oasys_util import read_surface_file

from syned.beamline.optical_elements.absorbers.filter import Filter
from wofrysrw.beamline.optical_elements.absorbers.srw_transmission import SRWTransmission

class SRWFilter(Filter, SRWTransmission):
    def __init__(self,
                 name="Undefined",
                 material="Be",
                 thickness=1e-3,
                 x_range=[-1e-3, 1e-3],
                 y_range=[-1e-3, 1e-3],
                 n_points_x=100,
                 n_points_y=100,
                 energy=15000,
                 thickness_error_profile=None,
                 scaling_factor=1.0):
        Filter.__init__(self, name=name, material=material, thickness=thickness)

        thickness_profile = numpy.ones((n_points_x, n_points_y))*self.get_thickness()

        if not thickness_error_profile is None:
            xx, yy, zz = read_surface_file(thickness_error_profile)
            interpolator = RectBivariateSpline(xx, yy, (zz.T)*scaling_factor)

            xx_t = numpy.linspace(x_range[0], x_range[1], n_points_x)
            yy_t = numpy.linspace(y_range[0], y_range[1], n_points_y)

            for i in range(n_points_x):
                for j in range(n_points_y):
                    thickness_profile[i, j] += interpolator.ev(xx_t[i], yy_t[j])

        attenuation_length, delta = get_absorption_parameters(material, energy)

        SRWTransmission.__init__(self,
                                 x_range=x_range,
                                 y_range=y_range,
                                 transmission_amplitudes=get_transmission_amplitudes(thickness_profile, attenuation_length),
                                 transmission_optical_path_difference=get_transmission_optical_path_difference(thickness_profile, delta),
                                 energy=energy)

def get_transmission_amplitudes(thickness_profile, attenuation_length):
    return numpy.exp(-0.5*thickness_profile/attenuation_length)

def get_transmission_optical_path_difference(thickness_profile, delta):
    return -delta*thickness_profile

def get_absorption_parameters(material, energy):
    energy_in_KeV = energy / 1000

    mu    = xraylib.CS_Total_CP(material, energy_in_KeV) # energy in KeV
    rho   = get_material_density(material)
    delta = 1 - xraylib.Refractive_Index_Re(material, energy_in_KeV, rho)

    return 0.01/(mu*rho), delta

def get_material_density(material_name):
    if material_name is None: return 0.0
    if str(material_name.strip()) == "": return 0.0

    try:
        compoundData = xraylib.CompoundParser(material_name)
        n_elements = compoundData["nElements"]
        if  n_elements == 1:
            return xraylib.ElementDensity(compoundData["Elements"][0])
        else:
            density = 0.0
            mass_fractions = compoundData["massFractions"]
            elements = compoundData["Elements"]
            for i in range (n_elements): density += xraylib.ElementDensity(elements[i])*mass_fractions[i]
            return density
    except:
        return 0.0
