import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq

# ==================================================
# == Setting up parameters ==

# Conversion Factors
IN = 0.0254                             # metres per inch
PSI = 6894.757293168                    # pascals per psi

# Tank Parameters
TANK_DENSITY = 2700                 # [kg/m^3]
TANK_LENGTH = 27 * IN               # [m]
TANK_OD = 8 * IN                    # [m]
TANK_THICKNESS = 0.125 * IN         # [m]
ULLAGE = 0.05
EMISSIVITY = 0.9                    # 0 is perfect mirror, 1 is perfect absorber
K_WALL = 100                        # 6061-T6 near LOX temp [W/m*K] (~155 at room temp)

# Atmosphere Parameters
P_ATMO = 101325                     # [Pa]
T_ATMO = 308                        # [K]
ATMOSPHERE = "Air"

# Fluid Parameters
TANK_PRESSURE = 490 * PSI           # [Pa]
T_LOX = 98                          # [K]
FLUID = "Oxygen"
# latent heat of vaporization (energy needed to turn 1kg of saturated liquid to vapor) [J/kg]
H_VAP = PropsSI("H", "P", P_ATMO, "Q", 1, FLUID) - PropsSI("H", "P", P_ATMO, "Q", 0, FLUID)

# Constants
G = 9.81                            # [m/s^2]
SIGMA = 5.67e-8                     # Stefan-Boltzmann constant [W/m^2*K^4]
# ==================================================

# ==================================================
# Geometry

def geometry():
    outer_radii = TANK_OD / 2   
    inner_radii = outer_radii - TANK_THICKNESS                                                              # [m]
    total_eff_height = TANK_LENGTH + 2 * outer_radii                                                        # including endcap [m]
    outer_surface_area = (2 * np.pi * outer_radii * TANK_LENGTH) + (4 * np.pi * outer_radii ** 2)                         # [m^2]
    internal_volume = (np.pi * (inner_radii ** 2) * TANK_LENGTH) + ((4/3) * np.pi * (inner_radii ** 3))     # [m^3]
    lox_mass = PropsSI("D", "T", T_LOX, "P", TANK_PRESSURE, FLUID) * (1-ULLAGE) * internal_volume
    internal_surface_area = (2 * np.pi * inner_radii * TANK_LENGTH) + (4 * np.pi * inner_radii ** 2)
    wall_mass = TANK_DENSITY * ((np.pi * (outer_radii ** 2 - inner_radii ** 2) * TANK_LENGTH) + ((4/3) * np.pi * (outer_radii ** 3 - inner_radii ** 3)))
    return total_eff_height, outer_surface_area, internal_volume, lox_mass, internal_surface_area, wall_mass

# ==================================================

# ==================================================
# Heat Transfer

# approximates the ratio of the buoyancy force to the viscous force acting on a fluid
def grashof_number(beta, delta_t, length, kinematic_viscosity):                                       # how much air's density change per kelvin (1/T for ideal gas) [1/K]
    return (G * beta * delta_t * length ** 3) / kinematic_viscosity ** 2

# ratio of momentum diffusion and heat diffusion
def prandtl_number(kinematic_viscosity, thermal_diffusivity):
    return kinematic_viscosity / thermal_diffusivity

# Rayleigh number is the ratio of buoyancy pushing the cold air down to the viscosity and heat diffusion resisting that motion
# Tells us the flow regime of the atmospheric air around, the circulation is natural convection and is the main way heat gets into the tank
def rayleigh_number(Gr, Pr):
    return Gr * Pr

# how many times more heat moves by convection than would move by pure conduction through still air over the same length
def nusselt_number(Pr, Ra):
    # use the Churchill Chu Correlation for a vertical cylinder
    # works best for a wide range of Ra values (10e-1 < Ra < 10e12)
    return (0.825 + ((0.387 * Ra ** (1/6)) / (1 + ((0.492 / Pr) ** (9/16))) ** (8/27))) ** 2

# convective heat transfer coefficient
def conv_heat_transfer_coeff(surface_t, total_eff_height):
    # convective heat transfer term
    film_temp = (surface_t + T_ATMO) / 2                                            # averaged temperature near the tank [K]
    density = PropsSI("D", "T", film_temp, "P", P_ATMO, ATMOSPHERE)                 # [kg/m^3]
    dynamic_viscosity = PropsSI("V", "T", film_temp, "P", P_ATMO, ATMOSPHERE)       # [Pa*s]
    kinematic_viscosity = dynamic_viscosity / density                           # [m^2/s]
    k = PropsSI("L", "T", film_temp, "P", P_ATMO, ATMOSPHERE)                       # thermal conducivity [W/m*K]
    cp_atmo = PropsSI("Cpmass", "T", film_temp, "P", P_ATMO, ATMOSPHERE)                 # specific heat [J/kg*K]
    thermal_diffusivity = k / (density * cp_atmo)                                    # how fast heat fissues through air [m^2/s]

    Gr = grashof_number(1 / film_temp, abs(T_ATMO - surface_t), total_eff_height, kinematic_viscosity)
    Pr = prandtl_number(kinematic_viscosity, thermal_diffusivity)
    Ra = rayleigh_number(Gr, Pr)
    Nu = nusselt_number(Pr, Ra)

    return (Nu * k) / total_eff_height

# radiation heat transfer coefficient
def rad_heat_transfer_coeff(surface_t):
    return EMISSIVITY * SIGMA * (surface_t ** 2 + T_ATMO ** 2) * (surface_t + T_ATMO)

def heat_leak(h_conv, h_rad, surface_t, outer_surface_area):
    Qdot = (h_conv + h_rad) * outer_surface_area * (T_ATMO - surface_t)
    return Qdot

# ==================================================

# ==================================================
# Boiloff Edge Case

def boil_off(Qdot):
    return Qdot / H_VAP

# ==================================================