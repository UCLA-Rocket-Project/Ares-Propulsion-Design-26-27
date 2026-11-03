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