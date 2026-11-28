import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq
from heat_transfer_functions import *

# ==================================================
# Heat Leak

# Resistance Network
def thermal_resistances(ext_h_conv, int_h_conv, h_rad, outer_surface_area, wetted_area):
    R_ext = 1 / ((ext_h_conv + h_rad) * outer_surface_area)
    R_in = 1 / (int_h_conv * wetted_area)
    return R_ext, R_in

# Rate of heat transfer (W)
# series path: ambient -> R_ext -> outer wall -> R_wall -> inner wall -> R_in -> LOX
def heat_leak(R_ext, R_wall, R_in):
    return (T_ATMO - T_LOX) / (R_ext + R_wall + R_in)

# ==================================================

# ==================================================
# Boiloff

# Rate of mass boiloff (kg/s)
def boil_off(Qdot):
    return Qdot / H_VAP

# ==================================================

# ==================================================
# Main function

def main():
    total_eff_height, outer_surface_area, internal_volume, lox_mass, wetted_area, liquid_height, wall_mass = geometry()
    print(f"total height: {total_eff_height:.3g} m")
    print(f"outer surface area: {outer_surface_area:.3g} m^2")
    print(f"internal volume: {internal_volume:.3g} m^3")
    print(f"lox mass: {lox_mass:.3g} kg")
    print(f"wetted area: {wetted_area:.3g} m^2")
    print(f"liquid height: {liquid_height:.3g} m")

    T_wall = T_LOX + WALL_DELTA_T
    print(f"assumed wall temperature: {T_wall:.1f} K")
    ext_h_conv = ext_conv_heat_transfer_coeff(T_wall, total_eff_height)
    h_rad = rad_heat_transfer_coeff(T_wall)
    int_h_conv = int_conv_heat_transfer_coeff(T_wall, liquid_height)
    print(f"convective heat transfer coefficient: {ext_h_conv:.3g} W/m^2*K")
    print(f"radiation heat transfer coefficient: {h_rad:.3g} W/m^2*K")
    print(f"internal convective heat transfer coefficient: {int_h_conv:.3g} W/m^2*K")

    R_ext, R_in = thermal_resistances(ext_h_conv, int_h_conv, h_rad, outer_surface_area, wetted_area)
    R_wall = wall_resistance()
    print(f"External thermal resistance: {R_ext:#.3g} K/W")
    print(f"Internal thermal resistance: {R_in:#.3g} K/W")
    print(f"Wall thermal resistance: {R_wall:#.3g} K/W")

    Qdot_lox = heat_leak(R_ext, R_wall, R_in)
    print(f"total heat leak into LOX: {Qdot_lox:.0f} W")

    mdot_boil = boil_off(Qdot_lox)
    print("\n=== Edge Case: Energy only towards phase change (boiloff) ===")
    print(f"LOX boil off rate: {mdot_boil:.3g} kg/s")

if __name__ == "__main__":
    main()

# ==================================================