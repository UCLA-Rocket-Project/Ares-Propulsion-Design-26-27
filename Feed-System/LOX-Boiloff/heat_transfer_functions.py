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
WALL_DELTA_T = 5                    # Empirically Determined [K]
# NIST cryogenic material properties for Aluminum 6061-T6 valid below 300 K to find cp
NIST_6061_CP = [46.6467, -314.292, 866.662, -1298.3, 1162.27, -637.795, 210.351, -38.3094, 2.96344]

# Atmosphere Parameters
P_ATMO = 101325                     # [Pa]
T_ATMO = 298                        # [K]
ATMOSPHERE = "Air"

# Fluid Parameters
TANK_PRESSURE = 490 * PSI           # [Pa]
T_LOX = 98                          # [K]
FLUID = "Oxygen"
# latent heat of vaporization (energy needed to turn 1kg of saturated liquid to vapor) [J/kg]
H_VAP = PropsSI("H", "P", P_ATMO, "Q", 1, FLUID) - PropsSI("H", "P", P_ATMO, "Q", 0, FLUID)
C_SF = 0.013                        # surface-fluid constant, estimation for LOX on Aluminum
N_ROHSENOW = 1.7                    # Prandtl exponent, 1.7 for anything but water

# Constants
G = 9.81                            # [m/s^2]
SIGMA = 5.67e-8                     # Stefan-Boltzmann constant [W/m^2*K^4]
# ==================================================

# ==================================================
# Tank Properties

# Important geometry values
def geometry():
    outer_radii = TANK_OD / 2   
    inner_radii = outer_radii - TANK_THICKNESS                                                              # [m]
    total_eff_height = TANK_LENGTH + 2 * outer_radii                                                        # including endcap [m]
    outer_surface_area = (2 * np.pi * outer_radii * TANK_LENGTH) + (4 * np.pi * outer_radii ** 2)                         # [m^2]
    internal_volume = (np.pi * (inner_radii ** 2) * TANK_LENGTH) + ((4/3) * np.pi * (inner_radii ** 3))     # [m^3]
    lox_mass = PropsSI("D", "T", T_LOX, "P", TANK_PRESSURE, FLUID) * (1-ULLAGE) * internal_volume
    internal_surface_area = (2 * np.pi * inner_radii * TANK_LENGTH) + (4 * np.pi * inner_radii ** 2)
    wall_mass = TANK_DENSITY * ((np.pi * (outer_radii ** 2 - inner_radii ** 2) * TANK_LENGTH) + ((4/3) * np.pi * (outer_radii ** 3 - inner_radii ** 3)))
    wetted_area = (1 - ULLAGE) * internal_surface_area
    liquid_height = (1 - ULLAGE) * total_eff_height
    return total_eff_height, outer_surface_area, internal_volume, lox_mass, wetted_area, liquid_height, wall_mass

# Tank Wall Thermal Resistance
def wall_resistance():
    R_cylinder = np.log((TANK_OD / 2) / (TANK_OD / 2 - TANK_THICKNESS)) / (2 * np.pi * K_WALL * TANK_LENGTH)
    R_endcap = (1 / (TANK_OD / 2 - TANK_THICKNESS) - 1 / (TANK_OD / 2)) / (4 * np.pi * K_WALL)
    return (1/R_cylinder + 1/R_endcap) ** -1

# Tank cp given temperature
def cp_wall(temperature):
    temperature = min(temperature, 300)
    log_T = np.log10(temperature)
    return 10 ** sum(a * log_T ** i for i, a in enumerate(NIST_6061_CP))

# ==================================================

# ==================================================
# Saturation Properties

# Find properties at saturated state, "Q", 0 -> saturated liquid, "Q", 1 -> saturated vapor
def sat_prop(P_fill):
    T_sat = PropsSI("T", "P", P_fill, "Q", 0, FLUID)                  # saturation temperature [K]
    rho_l = PropsSI("D", "P", P_fill, "Q", 0, FLUID)                  # saturated liquid density [kg/m^3]
    rho_v = PropsSI("D", "P", P_fill, "Q", 1, FLUID)                  # saturated vapor density [kg/m^3]
    h_fg = PropsSI("H", "P", P_fill, "Q", 1, FLUID) - PropsSI("H", "P", P_fill, "Q", 0, FLUID)   # latent heat [J/kg]
    sigma = PropsSI("I", "P", P_fill, "Q", 0, FLUID)                  # surface tension [N/m]
    mu_l = PropsSI("V", "P", P_fill, "Q", 0, FLUID)                   # liquid dynamic viscosity [Pa*s]
    cp_l = PropsSI("Cpmass", "P", P_fill, "Q", 0, FLUID)              # liquid specific heat [J/kg*K]
    Pr_l = PropsSI("Prandtl", "P", P_fill, "Q", 0, FLUID)             # liquid Prandtl number [-]

    return {
        "T_sat": T_sat,
        "rho_l": rho_l,
        "rho_v": rho_v,
        "h_fg": h_fg,
        "sigma": sigma,
        "mu_l": mu_l,
        "cp_l": cp_l,
        "Pr_l": Pr_l,
    }

# ==================================================

# ==================================================
# Film Properties

# Fluid properties at the film temperature between a wall and a reference temperature
def film_prop(wall_t, target_t, pressure, fluid):
    film_temp = (wall_t + target_t) / 2                                         # averaged temperature in the film [K]
    density = PropsSI("D", "T", film_temp, "P", pressure, fluid)                # [kg/m^3]
    dynamic_viscosity = PropsSI("V", "T", film_temp, "P", pressure, fluid)      # [Pa*s]
    k = PropsSI("L", "T", film_temp, "P", pressure, fluid)                      # thermal conductivity [W/m*K]
    cp = PropsSI("Cpmass", "T", film_temp, "P", pressure, fluid)                # specific heat [J/kg*K]
    return film_temp, density, dynamic_viscosity, k, cp

# ==================================================

# ==================================================
# Convective Heat Transfer

# approximates the ratio of the buoyancy force to the viscous force acting on a fluid
def grashof_number(beta, delta_t, length, kinematic_viscosity):
    return (G * beta * delta_t * length ** 3) / kinematic_viscosity ** 2

# ratio of momentum diffusion and heat diffusion
def prandtl_number(kinematic_viscosity, thermal_diffusivity):
    return kinematic_viscosity / thermal_diffusivity

# Rayleigh number is the ratio of buoyancy driving the fluid along the wall to the viscosity and heat diffusion resisting that motion
# Tells us the natural convection flow regime, used for both the ambient air outside the tank and the LOX inside the tank
def rayleigh_number(Gr, Pr):
    return Gr * Pr

# how many times more heat moves by convection than would move by pure conduction through still air over the same length
def nusselt_number(Pr, Ra):
    # use the Churchill Chu Correlation for a vertical cylinder
    # works best for a wide range of Ra values (1e-1 < Ra < 1e12)
    return (0.825 + ((0.387 * Ra ** (1/6)) / (1 + ((0.492 / Pr) ** (9/16))) ** (8/27))) ** 2

# convective heat transfer coefficient (air <-> outer wall)
def ext_conv_heat_transfer_coeff(wall_t, total_eff_height):
    # convective heat transfer term
    film_temp, density, dynamic_viscosity, k, cp_atmo = film_prop(wall_t, T_ATMO, P_ATMO, ATMOSPHERE)
    kinematic_viscosity = dynamic_viscosity / density                           # [m^2/s]
    thermal_diffusivity = k / (density * cp_atmo)                                    # how fast heat fissues through air [m^2/s]

    Gr = grashof_number(1 / film_temp, abs(T_ATMO - wall_t), total_eff_height, kinematic_viscosity)
    Pr = prandtl_number(kinematic_viscosity, thermal_diffusivity)
    Ra = rayleigh_number(Gr, Pr)
    Nu = nusselt_number(Pr, Ra)

    return (Nu * k) / total_eff_height

# convective heat transfer coefficient (inner wall <-> LOX)
def int_conv_heat_transfer_coeff(wall_t, liquid_height):
    film_temp, density, dynamic_viscosity, k, cp_lox = film_prop(wall_t, T_LOX, TANK_PRESSURE, FLUID)
    kinematic_viscosity = dynamic_viscosity / density
    thermal_diffusivity = k / (density * cp_lox)
    expansion_coefficient = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T", film_temp, "P", TANK_PRESSURE, FLUID)

    Gr = grashof_number(expansion_coefficient, abs(wall_t - T_LOX), liquid_height, kinematic_viscosity)
    Pr = prandtl_number(kinematic_viscosity, thermal_diffusivity)
    Ra = rayleigh_number(Gr, Pr)
    Nu = nusselt_number(Pr, Ra)
    return (Nu * k) / liquid_height

# ==================================================

# ==================================================
# Radiation Heat Transfer

def rad_heat_transfer_coeff(wall_t):
    return EMISSIVITY * SIGMA * (wall_t ** 2 + T_ATMO ** 2) * (wall_t + T_ATMO)

# ==================================================

