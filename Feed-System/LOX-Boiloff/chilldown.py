"""
Assumptions:
1. Tank only interface: Only the tank wall's thermal mass is considered; transfer lines,
   valves, and other connected hardware are excluded.
2. Vented fill at 1 atm: The tank is vented open to atmosphere during fill, so tank
   pressure is 1 atm and LOX is saturated at ~90.2 K (h_fg ~ 213 kJ/kg). With no
   subcooling, no heat can be absorbed as sensible heating of the liquid, making this
   the conservative boiloff case.
3. All heat reaching the liquid becomes boiloff: Vapor leaves through the vent at
   saturation temperature rather than warming up on the way out, which maximizes boiled
   mass (most conservative). Since the tank is open, all boiled GOX is lost.
4. Lumped wall: The tank wall is one uniform temperature with no gradient through or
   along the wall (valid in film boiling, approximate in the short nucleate phase).
5. Wetted tank: Inner wall is wetted up to the ullage line from t = 0; fill dynamics are
   ignored, but the full wall mass is still cooled.
6. Boiling regimes: Pool boiling progresses Film -> Transition -> Nucleate, with switch
   points calculated from correlations at 1 atm.
7. Heat leak temperature evaluation: Natural convection and radiation are evaluated at the
   current wall temperature, with air properties at the film temperature (Tw + T_atm)/2.
8. Still ambient air: Natural convection only, constant 298 K ambient; wind and frost are
   not modeled.
9. Chilldown end: Wall is considered chilled when it is within 5 K of LOX saturation
   temperature (~95 K).
"""

import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from heat_transfer_functions import *

# ==================================================
# Chilldown Boiling Regimes

# Minimum heat flux is when too little vapor is generated to hold the film up and thus where film boiling ends
def q_min(sat):
    # Berenson's Equation for Minimum Heat Flux [W/m^2]
   rho_l, rho_v = sat["rho_l"], sat["rho_v"]
   return 0.09 * rho_v * sat["h_fg"] * ((sat["sigma"] * G * (rho_l - rho_v)) / (rho_l + rho_v) ** 2) ** 0.25

# Maximum heat flux is when the superheat creates more bubbles and the vapor starts blocking incoming liquid, creating a film and is the peak of the boiling curve
def q_max(sat):
   # Zuber's Equation for Maximum Heat Flux [W/m^2]
   rho_l, rho_v = sat["rho_l"], sat["rho_v"]
   return 0.149 * sat["h_fg"] * rho_v * ((sat["sigma"] * G * (rho_l - rho_v)) / rho_v ** 2) ** 0.25

# Film boiling heat flux
# When the wall is too hot so heat crosses through a vapor layer which conducts poorly
def q_film(dT, sat):
   T_sat = sat["T_sat"]

   # Find vapor properties at halfway between tank wall and saturation temperature
   film_temp, density, dynamic_viscosity, k, cp = film_prop(T_sat + dT, T_sat, P_ATMO, FLUID)
   
   # Latent heat + vapor superheat [J/kg]
   h_fg_vap = sat["h_fg"] + 0.4 * cp * dT

   # Natural vapor bubble spacing on a large surface [m]
   capillary_length = np.sqrt(sat["sigma"] / (G * (sat["rho_l"] - sat["rho_v"])))

   # Film boiling heat transfer coefficient
   # Using Berenson's model for a large horizontal surface, predicting heat transfer near minimum heat flux point
   h_film = 0.425 * ((k ** 3 * density * (sat["rho_l"] - density) * G * h_fg_vap) / (dynamic_viscosity * dT * capillary_length)) ** 0.25

   return h_film * dT

# Rohsenow nucleate boiling correlation [W/m^2]
# Finds the heat flux during nucleate boiling without cap in order to find the Critical Heat Flux point later on
def q_rohsenow(dT, sat):
   dT = max(dT, 0)
   return (sat["mu_l"] * sat["h_fg"] * np.sqrt(G * (sat["rho_l"] - sat["rho_v"]) / sat["sigma"]) 
           * ((sat["cp_l"] * dT) / (C_SF * sat["h_fg"] * sat["Pr_l"] ** N_ROHSENOW)) ** 3)

# Nucleate boiling 
# When the liquid wets the wall and bubbles form, during this time the bubbles pulls cold liquid behind it and the heat flux rises
def q_nucleate(dT, sat):
   # Heat flux is capped given Zuber's Equation
   return min(q_rohsenow(dT, sat), q_max(sat))

# Leidenfrost Point, temperature where film boiling drops to the minimum heat flux
def dT_leidenfrost(sat, dT_upper):
   q_min_val = q_min(sat)
   return brentq(lambda dT: q_film(dT, sat) - q_min_val, 1, dT_upper)

# Critical Heat Flux Point, temperature where nucleate boiling reaches the CHF
def dT_CHF(sat, dT_upper):
   q_max_val = q_max(sat)
   return brentq(lambda dT: q_rohsenow(dT, sat) - q_max_val, 0.1, dT_upper)

# Transition boiling
# When the wall is partly wetted and partly covered in film, no reliable correlation so we use a log-linear interpolation between the two points
def q_transition(dT, boiling_constants):
   frac = (dT - boiling_constants["dT_chf"]) / (boiling_constants["dT_leid"] - boiling_constants["dT_chf"])          # 0 at CHF point, 1 at Leidenfrost point
   log_q = np.log10(boiling_constants["q_max"]) + frac * (np.log10(boiling_constants["q_min"]) - np.log10(boiling_constants["q_max"]))
   return 10 ** log_q

# Boiling regimes
def q_boil(dT, sat, boiling_constants):
   # temperature is greater than leidenfrost point -> film boiling
   # temperature is less than critical heat flux point -> nucleate boiling
   # temperature is in between -> transition boiling
   if dT > boiling_constants["dT_leid"]:
      return q_film(dT, sat)
   elif dT > boiling_constants["dT_chf"]:
      return q_transition(dT, boiling_constants)
   else:
      return q_nucleate(dT, sat)

# ==================================================

# ==================================================
# Tank Heat Transfer

# Heat Leak
# Evaluated at each wall temperature, assuming the tank is one lumped thermal mass
def heat_leak(T_w, total_eff_height, outer_surface_area):
   h_conv = ext_conv_heat_transfer_coeff(T_w, total_eff_height)      # [W/m^2*K]
   h_rad = rad_heat_transfer_coeff(T_w)                              # [W/m^2*K]
   return (h_conv + h_rad) * outer_surface_area * (T_ATMO - T_w)

# Chilldown
def chilldown(time, state, sat, boiling_constants, total_eff_height, outer_surface_area, wetted_area, wall_mass):
   # time and state parameters required for ODE solve_ivp
   T_w = state[0]
   dT = T_w - sat["T_sat"]                                     # Difference in temperature between current tank wall and saturation [K]

   # Heat Transfers
   Q_ext = heat_leak(T_w, total_eff_height, outer_surface_area)         # Heat into tank from ambient environment [W]
   Q_boil = q_boil(dT, sat, boiling_constants) * wetted_area            # Heat out from tank to LOX [W]

   # Wall cooling rate [K/s]
   # Rate of tank cooling down = total heat transfer / (wall mass * wall specific heat)
   dTw_dt = (Q_ext - Q_boil) / (wall_mass * cp_wall(T_w))

   # Boil off rate [kg/s]
   # Rate of liquid phase change = Heat from tank to LOX / LOX latent heat
   dm_dt = Q_boil / sat["h_fg"]

   return dTw_dt, dm_dt

# ==================================================

# ==================================================
# ODE Solver Events

# return value gives the distance to event and is checked at every timestep, the moment the sign flips to a negative the solver saves that timestep and the state it has
# .terminal tells the solver whether or not to stop if the event has been reached, in this case it stops at chilldown
# .direction tells the solver when to trigger the event, in this case the event is triggered when the return value returns negative


# Wall reaches the Leidenfrost point: film boiling ends
def leidenfrost_event(time, state, sat, boiling_constants, *args):
    return state[0] - (sat["T_sat"] + boiling_constants["dT_leid"])
leidenfrost_event.terminal = False
leidenfrost_event.direction = -1

# Wall reaches the CHF point: nucleate boiling begins
def chf_event(time, state, sat, boiling_constants, *args):
    return state[0] - (sat["T_sat"] + boiling_constants["dT_chf"])
chf_event.terminal = False
chf_event.direction = -1

# Wall is chilled: stop integrating
def chilled_event(time, state, sat, boiling_constants, *args):
    chilled_tolerance = 5
    return state[0] - (sat["T_sat"] + chilled_tolerance)
chilled_event.terminal = True
chilled_event.direction = -1

# ==================================================

# ==================================================
# Main Function

def main():
   # Find Boiling Regimes
   sat = sat_prop(P_ATMO)
   dT_initial = T_ATMO - sat["T_sat"]
   dT_leid = dT_leidenfrost(sat, dT_initial)
   boiling_constants = {
      "q_min": q_min(sat),
      "q_max": q_max(sat),
      "dT_leid": dT_leid,
      "dT_chf": dT_CHF(sat, dT_leid)
   }
   print("===== Boiling Regimes =====")
   print(f"Saturation Temperature: {sat["T_sat"]:.3g} K")
   print(f"q_min: {boiling_constants['q_min']:.3g} W/m^2")
   print(f"q_max: {boiling_constants['q_max']:.3g} W/m^2")
   print(f"Leidenfrost superheat: {boiling_constants['dT_leid']:.1f} K (wall {sat['T_sat'] + boiling_constants['dT_leid']:.1f} K)")
   print(f"CHF superheat: {boiling_constants['dT_chf']:.1f} K (wall {sat['T_sat'] + boiling_constants['dT_chf']:.1f} K)")
   print("==================================================\n")

   # Solving chilldown ODE
   total_eff_height, outer_surface_area, internal_volume, lox_mass, wetted_area, liquid_height, wall_mass = geometry()
   sol = solve_ivp(
      chilldown,                          # ODE function
      (0, 600),                           # time span [s], upper bound is conservative
      [T_ATMO, 0.0],                      # initial state, [tank wall = T_atmo, boiloff mass = 0]
      method="LSODA",                     # ODE solver method
      args=(sat, boiling_constants, total_eff_height, outer_surface_area, wetted_area, wall_mass),     # arguments to pass into ODE
      events=[leidenfrost_event, chf_event, chilled_event],       # specific check points to get data from
      max_step=1.0,                       # timestep [s]
      rtol=1e-6,                          # error tolerances
      atol=[1e-3, 1e-6]
   )
   if sol.status != 1:
      print(f"WARNING: chilldown criterion not reached: {sol.message}")
   t_leid = sol.t_events[0][0]         # [s]
   t_chf = sol.t_events[1][0]          # [s]
   t_chill = sol.t_events[2][0]        # [s]
   m_chill = sol.y[1, -1]              # [kg]
   print("===== Chilldown Time =====")
   print(f"Film boiling ends (Leidenfrost): {t_leid:.1f} s")
   print(f"Nucleate boiling begins (CHF): {t_chf:.1f} s")
   print(f"Chilldown time: {t_chill:.1f} s")
   print(f"Chilldown total mass loss: {m_chill:.2f} kg")
   print("==================================================\n")

   # Solving for boiloff Rate
   # Re-evaluate chilldown ODE at each saved time and take the dm/dt return value
   mdot = np.array(
      [chilldown(sol.t[i], sol.y[:, i], sat, boiling_constants, total_eff_height, outer_surface_area, wetted_area, wall_mass)[1] 
      for i in range(len(sol.t))])

   print("===== Boiloff Rate =====")
   print(f"Peak boiloff rate: {mdot[np.argmax(mdot)]:.3g} kg/s")
   print(f"Average boiloff rate: {(m_chill / t_chill):.3g} kg/s")
   print("==================================================\n")

   # Plot Boiloff Regime
   dT_range = np.linspace(1, T_ATMO - sat["T_sat"], 500)
   q_range = [q_boil(dT, sat, boiling_constants) for dT in dT_range]
   fig1, ax1 = plt.subplots()
   ax1.semilogy(sat["T_sat"] + dT_range, q_range)
   ax1.axvline(sat["T_sat"] + boiling_constants["dT_leid"], linestyle="--", label="Leidenfrost")
   ax1.axvline(sat["T_sat"] + boiling_constants["dT_chf"], linestyle=":", label="CHF")
   ax1.set_xlabel("Wall temperature [K]")
   ax1.set_ylabel("Heat flux [W/m^2]")
   ax1.set_title("Boiloff Regimes Heat Flux")
   ax1.legend()

   def mark_regimes(ax):
      ax.axvline(t_leid, linestyle="--", label=f"Film boiling ends ({t_leid:.1f} s)")
      ax.axvline(t_chf, linestyle=":", label=f"Nucleate boiling begins ({t_chf:.1f} s)")
      ax.axvline(t_chill, linestyle="-.", label=f"Chilldown complete ({t_chill:.1f} s)")

   # Plot Tank Chilldown Time
   fig2, ax2 = plt.subplots()
   ax2.plot(sol.t, sol.y[0])
   ax2.axhline(sat["T_sat"], color="gray", linestyle="-", linewidth=0.8, label=f"LOX T_sat ({sat['T_sat']:.1f} K)")
   mark_regimes(ax2)
   ax2.set_xlabel("Time [s]")
   ax2.set_ylabel("Wall temperature [K]")
   ax2.set_title("Tank Wall Chilldown Time")
   ax2.legend()

   # Plot Boiloff Rate
   fig3, ax3 = plt.subplots()
   ax3.plot(sol.t, mdot)
   mark_regimes(ax3)
   ax3.set_xlabel("Time [s]")
   ax3.set_ylabel("Boiloff Rate [kg/s]")
   ax3.set_title("Boiloff Rate")
   ax3.legend()

   plt.show()

if __name__ == "__main__":
   main()

# ==================================================
