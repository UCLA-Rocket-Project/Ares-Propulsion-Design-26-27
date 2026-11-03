import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq
from heat_transfer_functions import *

# ==================================================
# Thermal Expansion Edge Case

def thermal_expansion(Qdot_lox, internal_volume, lox_mass):
    cp_lox = PropsSI("Cpmass", "T", T_LOX, "P", TANK_PRESSURE, FLUID)
    Tdot = Qdot_lox / (lox_mass * cp_lox)

    # critical temperature, when the tank is filled with LOX
    def excess_volume(T):
        return lox_mass / PropsSI("D", "T", T, "P", TANK_PRESSURE, FLUID) - internal_volume
    
    T_crit = brentq(excess_volume, T_LOX, 150)
    time_T_crit = (T_crit - T_LOX) / Tdot
    return Tdot, T_crit, time_T_crit

# ==================================================

def main():
    total_eff_height, outer_surface_area, internal_volume, lox_mass, internal_surface_area, wall_mass = geometry()
    print(f"total height: {total_eff_height:.4f} m")
    print(f"outer surface area: {outer_surface_area:.4f}, m^2")
    print(f"internal volume: {internal_volume:.4f} m^3")
    print(f"lox mass: {lox_mass:.4f} kg")

    h_conv = conv_heat_transfer_coeff(T_LOX, total_eff_height)                     # [W/m^2*K]
    h_rad = rad_heat_transfer_coeff(T_LOX)                                       # [W/m^2*K]
    print(f"convective heat transfer coefficient: {h_conv:.4f} W/m^2*K")
    print(f"radiation heat transfer coefficient: {h_rad:.4f} W/m^2*K")
    Qdot_lox = heat_leak(h_conv, h_rad, T_LOX, outer_surface_area)                     # [W]
    print(f"total heat leak into LOX: {Qdot_lox} W")

    mdot_boil = boil_off(Qdot_lox)
    print("\n=== Edge Case 1: Energy only towards phase change (boiloff) ===")
    print(f"LOX boil off rate: {mdot_boil:.4f} kg/s")
    print("\n=== Edge Case 2: Energy only towards warming up LOX (thermal expansion overpressurization) ===")
    Tdot, T_crit, time_T_crit = thermal_expansion(Qdot_lox, internal_volume, lox_mass)
    print(f"Temperature change rate: {Tdot:.4f} K/s")
    print(f"Critical temperature: {T_crit:.4f} K")
    print(f"Time to critical temperature: {time_T_crit:.4f} s")
    

if __name__ == "__main__":
    main()