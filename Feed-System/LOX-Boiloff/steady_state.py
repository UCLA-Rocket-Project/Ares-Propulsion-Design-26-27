import numpy as np
import matplotlib.pyplot as plt
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
    print("===== Tank Geometry =====")
    print(f"total height: {total_eff_height:.3g} m")
    print(f"outer surface area: {outer_surface_area:.3g} m^2")
    print(f"internal volume: {internal_volume:.3g} m^3")
    print(f"lox mass: {lox_mass:.3g} kg")
    print(f"wetted area: {wetted_area:.3g} m^2")
    print(f"liquid height: {liquid_height:.3g} m")
    print("==================================================\n")

    T_wall = T_LOX + WALL_DELTA_T
    print("===== Heat Transfer Coefficients =====")
    print(f"saturated wall temperature: {T_wall:.1f} K")
    ext_h_conv = ext_conv_heat_transfer_coeff(T_wall, total_eff_height)
    h_rad = rad_heat_transfer_coeff(T_wall)
    int_h_conv = int_conv_heat_transfer_coeff(T_wall, liquid_height)
    print(f"convective heat transfer coefficient: {ext_h_conv:.3g} W/m^2*K")
    print(f"radiation heat transfer coefficient: {h_rad:.3g} W/m^2*K")
    print(f"internal convective heat transfer coefficient: {int_h_conv:.3g} W/m^2*K")
    print("==================================================\n")

    R_ext, R_in = thermal_resistances(ext_h_conv, int_h_conv, h_rad, outer_surface_area, wetted_area)
    R_wall = wall_resistance()
    print("===== Thermal Resistance Network =====")
    print(f"External thermal resistance: {R_ext:#.3g} K/W")
    print(f"Internal thermal resistance: {R_in:#.3g} K/W")
    print(f"Wall thermal resistance: {R_wall:#.3g} K/W")
    print("==================================================\n")

    Qdot_lox = heat_leak(R_ext, R_wall, R_in)
    mdot_boil = boil_off(Qdot_lox)
    print("===== Steady State Boiloff =====")
    print(f"total heat leak into LOX: {Qdot_lox:.0f} W")
    print(f"LOX boil off rate: {mdot_boil:.3g} kg/s")
    print("==================================================\n")

    # Temperature at each interface: walk the temperature drops outward from the LOX
    r_o = TANK_OD / 2                                       # outer wall radius [m]
    r_i = r_o - TANK_THICKNESS                              # inner wall radius [m]
    T_wi = T_LOX + Qdot_lox * R_in                          # inner wall surface [K]
    T_wo = T_wi + Qdot_lox * R_wall                         # outer wall surface [K]

    # Equivalent film thicknesses: a still layer of fluid with the same resistance, delta = k / h
    k_lox = film_prop(T_wi, T_LOX, TANK_PRESSURE, FLUID)[3]                 # LOX conductivity [W/m*K]
    k_air = film_prop(T_wo, T_ATMO, P_ATMO, ATMOSPHERE)[3]                  # air conductivity [W/m*K]
    delta_in = k_lox / int_h_conv                                           # inner LOX film [m]
    delta_out = k_air / (ext_h_conv + h_rad)                                # outer air film, conv + rad [m]

    # Radial temperature profile, region by region
    r_lox = np.linspace(0, r_i - delta_in, 50)
    T_lox = np.full_like(r_lox, T_LOX)                                      # bulk LOX, uniform
    r_film_in = np.linspace(r_i - delta_in, r_i, 50)
    T_film_in = T_LOX + (T_wi - T_LOX) * (r_film_in - (r_i - delta_in)) / delta_in      # linear across film
    r_wall = np.linspace(r_i, r_o, 50)
    T_wall_profile = T_wi + (T_wo - T_wi) * np.log(r_wall / r_i) / np.log(r_o / r_i)  # log across cylinder wall
    r_film_out = np.linspace(r_o, r_o + delta_out, 50)
    T_film_out = T_wo + (T_ATMO - T_wo) * (r_film_out - r_o) / delta_out               # linear across film
    r_amb = np.linspace(r_o + delta_out, r_o + 4 * delta_out, 20)
    T_amb = np.full_like(r_amb, T_ATMO)                                     # ambient

    regions = [
        (r_lox, T_lox, "Bulk LOX", "tab:blue"),
        (r_film_in, T_film_in, "LOX film", "tab:cyan"),
        (r_wall, T_wall_profile, "Aluminum wall", "tab:gray"),
        (r_film_out, T_film_out, "Air film (conv + rad)", "tab:orange"),
        (r_amb, T_amb, "Ambient air", "tab:red"),
    ]

    # Plot 1: full radial profile
    fig1, ax1 = plt.subplots()
    for r, T, label, color in regions:
        ax1.plot(r * 1000, T, color=color, linewidth=2, label=label)
    ax1.set_xlabel("Distance from tank centerline [mm]")
    ax1.set_ylabel("Temperature [K]")
    ax1.set_title(f"Steady-state radial temperature profile (Q = {Qdot_lox:.0f} W)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    # Plot 3: thermal resistance
    names = ["Air film (R_ext)", "Aluminum wall (R_wall)", "LOX film (R_in)"]
    values = [R_ext, R_wall, R_in]
    R_total = sum(values)
    fig3, ax3 = plt.subplots()
    bars = ax3.barh(names, values, color=["tab:orange", "tab:gray", "tab:cyan"])
    ax3.set_xscale("log")
    for bar, R in zip(bars, values):
        ax3.text(R * 1.15, bar.get_y() + bar.get_height() / 2,
                 f"{R:.2e} K/W ({100 * R / R_total:.2f}%)", va="center", fontsize=8)
    ax3.set_xlim(min(values) / 3, max(values) * 30)
    ax3.set_xlabel("Thermal resistance [K/W] (log scale)")
    ax3.set_title(f"Thermal resistance network (total {R_total:.3g} K/W)")

    plt.show()

if __name__ == "__main__":
    main()

# ==================================================