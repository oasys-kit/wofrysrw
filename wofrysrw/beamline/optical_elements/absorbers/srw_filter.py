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
                 thickness_error_profile=None):
        Filter.__init__(self, name=name, material=material, thickness=thickness)

        transmission_optical_path_difference = numpy.zeros((n_points_x, n_points_y))
        thickness_profile = numpy.full((n_points_x, n_points_y), self.get_thickness())

        if not thickness_error_profile is None:
            xx, yy, zz = read_surface_file(thickness_error_profile)
            interpolator = RectBivariateSpline(xx, yy, zz.T)

            xx_t = numpy.linspace(x_range[0], x_range[1], n_points_x)
            yy_t = numpy.linspace(y_range[0], y_range[1], n_points_y)

            for i in range(n_points_x):
                for j in range(n_points_y):
                    thickness_error = interpolator.ev(xx_t[i], yy_t[i])
                    transmission_optical_path_difference[i, j] += thickness_error
                    thickness_profile += thickness_error

        transmission_amplitudes = getAmplitudeTransmittance(thickness_profile, alpha=getLinearAbsorptionCoefficient(self.get_material(), energy))

        SRWTransmission.__init__(self,
                                 x_range=x_range,
                                 y_range=y_range,
                                 transmission_amplitudes=transmission_amplitudes,
                                 transmission_optical_path_difference=transmission_optical_path_difference,
                                 energy=energy)

def getAmplitudeTransmittance(thickness, alpha):
    return numpy.sqrt(numpy.exp(-alpha * thickness * 100))

def getLinearAbsorptionCoefficient(chemical_formula, energy):
    mu = xraylib.CS_Total_CP(chemical_formula, energy/1000) # energy in KeV
    rho = getMaterialDensity(chemical_formula)

    return mu*rho

def getMaterialDensity(material_name):
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
