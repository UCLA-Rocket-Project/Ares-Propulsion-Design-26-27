function results = claude_boiloff(varargin)
%BOILOFF_V3  boiloff_v2 + PRESSURIZATION. Stratified (non-equilibrium) model of the vent / top-off fill of an uninsulated
%  cryogen tank from a dewar (LOX or LN2; 304 SS or 6061-T6 Al, set in section 1). Same tank, lines,
%  vent, wall and procedure as boiloff.m, but the
%  tank fluid is split into three nodes with their own temperatures:
%
%    * ULLAGE VAPOR   mass m_v, internal energy U_v -> T_v, P from the real gas (can be superheated)
%    * SURFACE LAYER  top dz_surf of liquid, mass m_s, temperature T_s
%    * BULK LIQUID    the rest of the liquid, mass m_b, temperature T_b (can be subcooled)
%
%  The tank pressure is the ullage pressure. The liquid-vapor interface sits at Tsat(P):
%    evaporation / condensation  m_ev = (Q_surface->interface + Q_vapor->interface)/h_fg
%      Q_surface->interface = h_li*A*(T_s - Tsat),  h_li = h_evap if T_s > Tsat else h_cond
%      Q_vapor->interface   = h_vi*A*(T_v - Tsat)
%  Surface layer <-> bulk: weak conduction when the warmer liquid is on top (stable), strong
%  mixing when the surface is colder than the bulk. Surface layer mass relaxes to its nominal
%  thickness (refilled from the bulk in tau_up, drains back to the bulk in tau_down).
%  Incoming liquid flashes to the tank pressure: the flash vapor rises through the bulk (collapsing
%  if the bulk is subcooled, like wall bubbles) and the rest reaches the ullage; the liquid
%  rises to the surface layer if it is warmer than the bulk, otherwise it joins the bulk.
%  Wall boiling heats the bulk while it is subcooled (bubbles collapse in it) and turns into
%  vapor once the bulk is within dT_sub_vap of saturation. A superheated bulk flashes (tau_boil).
%  Liquid leaving through the dip tube comes from the surface layer; backflow to the dewar
%  comes from the bulk.
%
%  PRESSURIZATION (new in v3, after the sealed hold until fire):
%    GN2 K bottles -> regulator (outlet droops with flow) -> one 1/32 in solenoid orifice -> manifold
%    -> check valve + branch -> LOX tank and fuel tank. The manifold pressure is where the solenoid flow
%    equals the two branch flows. Bottles and the fuel-tank ullage are real-gas nitrogen with wall heat.
%    In the LOX tank the ullage becomes an O2 + N2 mixture (Dalton, one temperature):
%      * the O2 evaporates / condenses at the temperature where its vapor pressure over the surface
%        liquid equals its PARTIAL pressure (Raoult), so warm N2 does not stop O2 boil-off;
%      * boiling at the wall, flashing and bubble survival use Tsat at the TOTAL pressure, so as the
%        pressure climbs the liquid becomes subcooled and wall heat only warms it (above ~590 psig
%        there is no boiling at all: 800 psig is above the O2 critical pressure, 717 psig);
%      * N2 dissolves / condenses into the surface layer (Raoult equilibrium, rate set by beta_abs),
%        releasing its condensation heat there; dissolved N2 is tracked in the surface and the bulk;
%      * the liquid is compressed by the pressure (CoolProp compressibility), freeing ullage volume.
%    Outputs the LOX-only liquid mass at fire, the dissolved N2, the vapor masses and the subcooling.
%
%  USAGE
%    boiloff_v3                              % run with the inputs in section 1
%    boiloff_v3('T_dewar_liq',100)           % any p.<name> input by name
%    boiloff_v3('T_wall0',linspace(130,260,20))   % starting wall profile, bottom -> top
%    boiloff_v3('cryogen','LN2','material','6061Al')   % LN2 into an aluminum tank
%    boiloff_v3('t_press',522,'P_reg_psig',816,'OF',1,'m_fuel',22.5,'meas_press_file','Data/press_1s_hf4.csv')
%
%  REQUIREMENTS: MATLAB R2021a+, Python with CoolProp visible to MATLAB (see boiloff.m).
%  Property tables are cached per cryogen (lox_tables_v3.mat, ln2_tables_v3.mat) next to this file
%  (first run for each cryogen ~5 min: v2's tables plus the pressurant tables).
%
%  KNOWN LIMITATIONS
%    1. One temperature for all the ullage vapor, so cold gas entering at the bottom (pre-chill)
%       cools the whole wall evenly. Use a T_wall0 profile if you expect the bottom to be colder.
%    2. Interface and mixing coefficients (h_vi, h_evap, h_cond, h_sb_*) are order-of-magnitude
%       values, not correlations. They control how fast the pressure follows the liquid.
%    3. The dewar liquid is saturated, and with dewar_sag its pressure falls as liquid is drawn.
%       V_dewar and dewar_fill are stand-ins that set how fast it falls.
%    4. Tank treated as a vertical cylinder; end caps neglected.
%    5. Pressurization: the N2 absorption rate (beta_abs) is an estimate; at 800 psig (above the O2
%       critical pressure) there is no true liquid-vapor interface, so the interface model is an
%       engineering approximation. Bottle size, fuel load and PT calibrations are assumptions (see inputs).

%% 1. USER INPUTS ----------------------------------------------------------
psi2Pa = 6894.757;  in2m = 0.0254;  ft2m = 0.3048;

% ---- Material and cryogen ----
% Everything that depends on them follows from these: fluid properties (CoolProp), wall
% conductivity and specific heat (NIST cryogenic fits), density and emissivity.
p.cryogen       = 'LOX';    % 'LOX' or 'LN2'
p.material      = '304SS';  % tank wall: '304SS' or '6061Al' (6061-T6)
p.line_material = '304SS';  % uninsulated fill-line tubes: '304SS', '6061Al', or 'tank' (same as the tank)

% ---- Tank ----
p.D        = 7.981*in2m;     % tank inner diameter                                 [m]
p.V_tank   = 1938.1*in2m^3;  % total internal tank volume                          [m^3]
p.t_w      = 0.322*in2m;     % tank wall thickness                                 [m]
p.V_dip    = 1370.4*in2m^3;  % liquid volume when the level reaches the dip-tube inlet [m^3]
p.N_wall   = 20;             % number of axial wall bands
p.k_w      = NaN;            % OVERRIDES wall conductivity for axial conduction [W/m/K]. NaN = NIST curve k(T) of
                             % p.material (304: 8.7 at 90 K, 15.3 at 300 K; 6061: 92 / 155); a number = constant

% ---- Vent (dip tube) ----
p.vtube_OD   = 0.25*in2m;    % vent / dip tube outer diameter (-4)                   [m]
p.vtube_wall = 0.049*in2m;   % vent tube wall as installed (0.152 in bore)           [m]
p.Cv_vent_valve = 9.7;       % vent valve Cv (same as the fill valve)
p.ID_vent_valve = NaN;       % vent valve bore [m]. NaN = no reducer losses; set it (e.g. 0.37*in2m)
                             % to add the expansion + contraction from the 0.152 in tube into the valve
p.L_vent     = 18.68*in2m;   % vent tube total length: 11.73 in inside the tank (dip tube) + 6.75 in
                             % above it. All of it carries friction. NaN = tank length above the dip inlet [m]

% ---- Fill path: dewar tube -> globe valve -> hose -> 1/2 in tube -> ball valve -> 1/2 in -> 3/4 in -> tank ----
p.L_hose   = 5*ft2m;         % insulated cryogenic flex hose (-8 AN), adiabatic      [m]
p.ID_hose  = 0.5*in2m;       % ASSUMED hose bore                                 [m]
p.eps_hose = 1.5e-6;         % hose roughness                                        [m]
p.OD1 = 0.5*in2m;  p.t1 = 0.065*in2m;  p.L1 = 1.5*ft2m;   % 1/2 in tube, dewar side of the ball valve
p.OD2 = 0.5*in2m;  p.t2 = 0.065*in2m;  p.L2 = 1.0*ft2m;   % 1/2 in tube, ball valve -> 3/4 in tube
p.OD3 = 0.75*in2m; p.t3 = 0.065*in2m;  p.L3 = 1.5*ft2m;   % 3/4 in tube -> tank bottom
p.eps_ln     = 1.5e-6;       % drawn tube roughness                                  [m]
p.K_ent      = 0.5;          % entrance loss (square-edged inlet)
p.Cv_valve   = 9.7;          % fill ball valve Cv, fully open
p.valve_open = 0.9;          % fraction open (0-1)
p.valve_curve = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1; ...
                 0 .02 .05 .09 .15 .23 .34 .48 .65 .83 1];   % typical full-port ball valve
p.Cv_dewar   = 1.1;          % dewar liquid globe valve (RegO 3/8 in: Cv 1.10). NaN = not modeled
p.L_dw  = 1.6;               % <-- ASSUMED dewar withdrawal tube + external tube length   [m]
p.ID_dw = 0.305*in2m;        % 3/8 x 0.035 in tube (real bore)                         [m]
p.K_fit      = 0;            % extra lumped loss on the 1/2 in tube velocity
p.backflow   = true;         % tank liquid pushed back to the dewar when P_tank > P_dewar

% ---- Operating conditions ----
p.P_d_psig = 130;            % dewar pressure at the START of the main fill [psig] (hf4). No pressure builder,
                             % so it sags as liquid is drawn; 120 sags to ~109 by the hf2 plateau.
p.T_dewar_liq = NaN;         % dewar liquid temperature [K]. NaN = saturated (hf2: no pressure builder).
p.T_amb    = 286%277;            % ambient temperature                                   [K]
p.P_amb    = 92389.75;       % ambient pressure                                      [Pa]
p.T_wall0  = NaN;            % initial wall temperature [K]: NaN = T_amb, a scalar, or N_wall values
                             % from bottom to top (applied before the pre-chill)
p.T_ln1_0  = NaN;            % tube 1 start temperature; NaN = dewar liquid temperature [K]

% ---- Procedure ----
p.n_fill   = 3;              % number of fill phases (hf4: fill + 3 refills)
p.t_vent   = 120;            % vent-only phase durations, scalar or one per vent (hf4: 120 s) [s]
p.tA_min   = 0%[440 15 20 15]; % minimum duration of each fill phase (hf4: T+3346-3786, then refills) [s]
p.tA_extra = 7;      % keep filling this long after liquid first reaches the vent (Assume 5 seconds for reaction time) [s]
p.final_vent = true;         % after the last fill, vent until liquid stops coming out
p.t_hold   = 50;            % then vent AND fill closed (sealed tank) until pressurization [s] (hf4 ~50 s)
                             % closed ~T+2622, pressure rose 67 -> 80 psig until pressurization ~T+2665

% ---- Matching Ablative 2 Dewar Swap ----
p.t_prechill  = 1045;           % optional pre-chill from a low-pressure source [s]; 0 = none (hf4).
                             % e.g. 1045 s from a 30 psig source holding the tank at 24 psig
p.P_src_psig  = 30;          % pre-chill source pressure                             [psig]
p.P_hold_psig = 24;          % tank pressure held during pre-chill                   [psig]
p.prechill_x  = 1;           % quality arriving during pre-chill (1 = gas)
p.tA_max   = 3600;           % safety cap on one fill phase                          [s]
p.tC_max   = 1200;           % safety cap on the final vent                          [s]

% ---- Pressurization (after the sealed hold, until fire) ----
% One low-flow solenoid feeds BOTH tanks from the GN2 K bottles through a regulator; each tank has a check
% valve on its pressurant branch. Cd_sol = 0.6 and reg_droop = 14.5 psi per g/s reproduce the measured
% hf4 / hf2 tank-pressure rise within ~10-15 psi RMS (pressurization figure). Defaults below are hf2.
p.t_press      = 543;          % solenoid open until fire [s] (hf2 543, hf4 522); 0 = no pressurization
p.press_gas    = 'Nitrogen';   % pressurant (CoolProp name): 'Nitrogen' or 'Helium'
p.P_bot0_psig  = 2200;         % K-bottle pressure when the solenoid opens                     [psig]
p.n_bot        = 2;            % K bottles in parallel
p.V_bot1       = 49.9e-3;      % <-- ASSUMED internal (water) volume of one K bottle            [m^3]
p.D_bot = 0.2286;  p.L_bot = 1.30;   % K-bottle OD and height (gas-to-wall heat transfer)     [m]
p.P_reg_psig   = 794;          % regulator outlet at no flow (data: hf2 794, hf4 816 psig)      [psig]
p.reg_droop    = 14.5;         % regulator outlet droop per unit flow (data: ~60 psi at first)  [psi per g/s]
p.d_sol        = in2m/32;      % solenoid orifice diameter                                       [m]
p.Cd_sol       = 0.6;          % solenoid discharge coefficient (fitted to the hf2 / hf4 rise)
p.d_branch     = 0.245*in2m;  p.Cd_branch = 0.6;   % each tank's branch (3/8 in tube)        [m]
p.P_crack_psig = 0;            % check-valve cracking pressure                                   [psig]
p.beta_abs     = 2e-5;         % N2 absorption into the LOX surface, liquid-side mass
                               %     transfer coefficient [m/s]. 1e-5 to 4e-5 all fit the hf4 bottle drop
p.beta_sb      = 1e-6;         % dissolved N2 surface -> bulk exchange [m/s] (N2-rich liquid is lighter)
% ---- Fuel tank (shares the pressurant manifold) ----
p.V_fuel_tank  = 1938.1*in2m^3;  p.D_fuel_tank = 7.981*in2m;   % same tank as the LOX tank
p.m_lox_target = 56.5*0.45359237;  % target LOX load (56.5 lb)                                  [kg]
p.OF           = 1.3;          % mixture ratio (hf4 used 1.0)
p.m_fuel       = NaN;          % fuel loaded [kg]; NaN = m_lox_target/OF (19.7 kg). The early fuel-tank
                               % rise in hf4 implies ~22-23 kg (OF 1 on the LOX actually loaded)
p.rho_fuel     = 860;          % 75 % ethanol / 25 % water by mass, ~13 C                       [kg/m^3]

% ---- Stratification (interface, surface layer, bulk) ----
% Magnitude matters not value (2x vs 5x makes little difference) but fitted
% to match fill HF data, should be changed for LN2
p.dz_surf     = 0.01;        % surface layer thickness                               [m]
p.h_vi        = 20;          % vapor -> interface heat transfer coefficient          [W/m2/K]
p.h_evap      = 2000;        % superheated surface layer -> interface (evaporation)  [W/m2/K]
p.h_cond      = 3000;        % interface -> subcooled surface layer (condensation)    [W/m2/K]
                             % set so the sealed-tank rise is ~0.3 psi/s as in hf2 (200 gives ~0.7)
p.h_sb_stable = 30;          % surface <-> bulk, warm liquid on top (stable)         [W/m2/K]
p.h_sb_mix    = 1000;        % surface <-> bulk, cold liquid on top (overturning)     [W/m2/K]
p.tau_up      = 5;           % surface layer refill from the bulk                    [s]
p.tau_down    = 300;         % excess surface liquid draining into the bulk          [s]
p.tau_boil    = 2;           % flashing time of a superheated bulk                   [s]
p.dT_sub_vap  = 1.5;         % bulk subcooling below which wall bubbles reach the ullage [K]

% ---- Bubbles in the liquid / dip tube ----
p.swell    = false;           % bubble (void fraction) model on/off. On: tank ends ~21.3 kg (hf2: 21.6
                             % kg expelled); off: ~23.7 kg = bubble-free liquid filling V_dip
p.C0       = 1.2;  p.Ugj_coef = 1.41;
p.dh_dip   = 0.005;          % height over which the dip-tube inlet goes from vapor to liquid [m]

% ---- Dewar ----
p.dewar_sag  = true;         % true = no pressure builder: P_dewar falls as liquid is drawn off (hf4: 108 -> 92 psig)
p.V_dewar    = 0.230;        % dewar total volume (230 L)                             [m^3]
p.dewar_fill = 0.75;         % liquid volume fraction at the start of the fill (just under half)
                             % NOTE: the sag model keeps the dewar's liquid and vapor at one temperature,
                             % which is the SLOWEST possible sag (~5-6 psi here vs ~16 psi seen in hf4).

% ---- Heat transfer constants & Overrides----
p.q_ext_mult = 1.0;  p.C_sf = 0.013;  p.dTmin_ovr = NaN;
p.rho_w    = NaN;    % tank wall density [kg/m3]. NaN = material default (304: 7900, 6061: 2700)
p.eps      = NaN;    % tank outside emissivity. NaN = material default (304: 0.3, 6061: 0.1; estimates,
                     % surface finish and frost change these a lot)
p.eps_line = NaN;    % fill-line outside emissivity. NaN = default of p.line_material
p.gam_v  = NaN;     % vent gas exponent. NaN = live from CoolProp: isentropic exponent of the ullage gas at
                    % (T_v, P), the value that sets its real sound speed (~1.39-1.40 here). A number = fixed.
p.f_vent = NaN;     % vent tube Darcy friction factor. NaN = live: Churchill f at the vent's own Reynolds
                    % number (G*D/mu, mu of the ullage gas), solved with the flow. A number = fixed.

% ---- Plots / measured overlay ----
p.make_plots  = true;
p.progress_dt = 10;          % print a progress line every this many seconds of real (wall-clock) time; 0 = off
p.show_meas   = false;
p.meas_file   = 'Data/ox_pt_1s_counts_hf2.csv';  
p.meas_zero   = 477.5;  p.meas_cal = 0.1907;  p.meas_t_fill = 1400;   % hf4: main fill starts at T+3346
p.meas_mass   = 22.85%22.3;      % from HF [kg]
p.meas_press_file = 'Data/press_1s_hf2.csv';   % measured pressurization (t from solenoid open, LOX, fuel,
                                               % bottle, regulator psig); '' = none. hf4: press_1s_hf4.csv

for i = 1:2:numel(varargin)
    if ~isfield(p, varargin{i}), error('boiloff_v3:input','Unknown input p.%s.', varargin{i}); end
    p.(varargin{i}) = varargin{i+1};
end

%% 2. DERIVED CONSTANTS AND PROPERTY TABLES --------------------------------
here = fileparts(mfilename('fullpath'));
p.cryo  = cryo_db(p.cryogen);   p.fluid = p.cryo.coolprop;
if strcmpi(p.line_material, 'tank'), p.line_material = p.material; end
p.mw  = mat_db(p.material);      p.mln = mat_db(p.line_material);
if isnan(p.rho_w),    p.rho_w    = p.mw.rho;  end
if isnan(p.eps),      p.eps      = p.mw.eps;  end
if isnan(p.eps_line), p.eps_line = p.mln.eps; end
fprintf('%s into a %s tank (%s fill lines)\n', p.cryo.name, p.mw.name, p.mln.name);
p.g = 9.80665;  p.sigSB = 5.670374e-8;
p.A_c   = pi*p.D^2/4;  p.L_t = p.V_tank/p.A_c;  p.h_dip = p.V_dip/p.A_c;
p.dz    = p.L_t/p.N_wall;  p.z0 = (0:p.N_wall-1).'*p.dz;  p.Ab = pi*p.D*p.dz;
p.mb    = p.rho_w*p.t_w*p.Ab;  p.C_ax = pi*p.D*p.t_w/p.dz;
% Axial conduction uses the conductivity integral theta(T) = int k dT, so Q = C_ax*(theta(T2) - theta(T1))
% is exact for k(T) even with a large temperature step between bands.
Tk = (1:0.5:400).';
if isnan(p.k_w), kk = k_mat(Tk, p.mw); else, kk = p.k_w*ones(size(Tk)); end
p.theta = griddedInterpolant(Tk, cumtrapz(Tk, kk), 'pchip', 'linear');
p.ID_vent = p.vtube_OD - 2*p.vtube_wall;  p.A_v = pi*p.ID_vent^2/4;
for s = '123'
    p.(['ID' s])   = p.(['OD' s]) - 2*p.(['t' s]);
    p.(['m_ln' s]) = p.mln.rho*pi/4*(p.(['OD' s])^2 - p.(['ID' s])^2)*p.(['L' s]);
end
p.A_in    = pi*p.ID1^2/4;
Cv_eff    = p.Cv_valve*interp1(p.valve_curve(1,:), p.valve_curve(2,:), p.valve_open);
p.K_valve = 894*(p.ID1/in2m)^4/Cv_eff^2;
p.P_d     = p.P_d_psig*psi2Pa + p.P_amb;

TAB_VERSION = 1;
cacheFile = fullfile(here, p.cryo.cache);  tabRaw = [];
if exist(cacheFile, 'file')
    S = load(cacheFile);
    if isfield(S,'tabRaw') && S.tabRaw.version == TAB_VERSION && strcmp(S.tabRaw.fluid, p.fluid) ...
            && isfield(S.tabRaw,'press') && strcmpi(S.tabRaw.press, p.press_gas)
        tabRaw = S.tabRaw;
    end
end
if isempty(tabRaw)
    fprintf('Building %s + %s property tables with CoolProp (one time, ~5 min)...\n', p.cryo.name, p.press_gas);
    tabRaw = build_tables(p.fluid, p.cryo, p.press_gas);  tabRaw.version = TAB_VERSION;  tabRaw.fluid = p.fluid;
    tabRaw.press = p.press_gas;
    save(cacheFile, 'tabRaw');
end
if ~isfield(tabRaw, 'gasKap')
    fprintf('Adding the gas isentropic exponent to the property cache (one time)...\n');
    tabRaw.gasKap = build_kappa(tabRaw, p.fluid);  save(cacheFile, 'tabRaw');
end
[p.tab, p.air] = make_interpolants(tabRaw);

p.Tsat_d = Tsat_of_P(p.P_d, p.tab);
p.sd     = sat_props(p.Tsat_d, p.tab);
if isnan(p.T_dewar_liq), p.T_dl = NaN;  p.s_dl = p.sd;
else, p.T_dl = min(p.T_dewar_liq, p.Tsat_d);  p.s_dl = sat_props(p.T_dl, p.tab);  p.dewar_sag = false; end
if isnan(p.T_ln1_0), p.T_ln1_0 = p.s_dl.T; end
p.P_hold = p.P_hold_psig*psi2Pa + p.P_amb;
s_src    = sat_props(Tsat_of_P(p.P_src_psig*psi2Pa + p.P_amb, p.tab), p.tab);
p.h_src  = s_src.h_l + p.prechill_x*(s_src.h_g - s_src.h_l);
p.m_s0   = p.A_c*p.dz_surf*sat_props(Tsat_of_P(p.P_amb, p.tab), p.tab).rho_l;   % nominal surface layer mass

p = setup_flows(p);

% ---- pressurant ----
p.same_gas = strcmpi(p.press_gas, p.fluid);              % e.g. LN2 pressurized with GN2: one species
p.press_dissolves = ~p.same_gas && strcmpi(p.press_gas, 'Nitrogen');   % helium: negligible solubility
p.M_c = cprops('M','T',300,'P',1e5,p.fluid);  p.M_p = cprops('M','T',300,'P',1e5,p.press_gas);   % [kg/mol]
p.R_p = 8.314462618/p.M_p;
cp0 = cprops('CP0MASS','T',p.T_amb,'P',1e5,p.press_gas);  p.gam_p = cp0/(cp0 - p.R_p);   % ideal-gas gamma
p.P_crit = cprops1('Pcrit', p.fluid);
p.P_bot0 = p.P_bot0_psig*psi2Pa + p.P_amb;  p.P_reg = p.P_reg_psig*psi2Pa + p.P_amb;
p.k_droop = p.reg_droop*psi2Pa*1000;                     % [Pa per kg/s]
p.CdA_sol = p.Cd_sol*pi*p.d_sol^2/4;  p.CdA_br = p.Cd_branch*pi*p.d_branch^2/4;
p.P_crack = p.P_crack_psig*psi2Pa;
p.V_bot = p.n_bot*p.V_bot1;  p.A_bot = p.n_bot*pi*p.D_bot*p.L_bot;
if isnan(p.m_fuel), p.m_fuel = p.m_lox_target/p.OF; end
p.D_fu = p.D_fuel_tank;  p.A_fuc = pi*p.D_fu^2/4;  p.V_fu = p.V_fuel_tank - p.m_fuel/p.rho_fuel;
if p.V_fu <= 0, error('boiloff_v3:fuel','%.1f kg of fuel does not fit in the fuel tank.', p.m_fuel); end
if p.t_press > 0
    fprintf('Pressurization: %s, %d x %.1f L bottles at %.0f psig, regulator %.0f psig, %.3f in orifice (Cd %.2f)\n', ...
        p.press_gas, p.n_bot, p.V_bot1*1e3, p.P_bot0_psig, p.P_reg_psig, p.d_sol/in2m, p.Cd_sol);
    fprintf('Fuel tank: %.1f kg fuel, ullage %.2f L\n', p.m_fuel, p.V_fu*1e3);
end
fprintf('Inflow: max %.3f kg/s (choked, full dewar pressure)\n', p.in_Gmax*p.A_in);

sdw = p.sd;
mld = p.dewar_fill*p.V_dewar*sdw.rho_l;  mgd = (1-p.dewar_fill)*p.V_dewar*sdw.rho_g;
Md0 = mld + mgd;  Ud0 = mld*sdw.u_l + mgd*sdw.u_g;
p.pdew = struct('tab',p.tab,'V_tank',p.V_dewar,'fluid',p.fluid);

%% 3. INITIAL CONDITIONS ---------------------------------------------------
rho0 = cprops('D','T',p.T_amb,'P',p.P_amb,p.fluid);
M0   = rho0*p.V_tank;  U0 = M0*cprops('U','T',p.T_amb,'P',p.P_amb,p.fluid);
T0l  = Tsat_of_P(p.P_amb, p.tab);
Tw0  = p.T_wall0;  if all(isnan(Tw0)), Tw0 = p.T_amb; end
Tw0  = Tw0(:).*ones(p.N_wall,1);
% Pressurant states: bottles full at P_bot0, T_amb; fuel ullage starts as ambient gas (air taken as N2)
rb = cprops('D','T',p.T_amb,'P',p.P_bot0,p.press_gas);  ub = cprops('U','T',p.T_amb,'P',p.P_bot0,p.press_gas);
rf = cprops('D','T',p.T_amb,'P',p.P_amb,p.press_gas);   uf = cprops('U','T',p.T_amb,'P',p.P_amb,p.press_gas);
% x = [m_v; U_v; m_s; T_s; m_b; T_b; T_ln1; T_ln2; T_ln3; m_vent; m_in; M_dew; U_dew; alpha; T_wall(1:N);
%      m_pv; m_ps; m_pb; m_bot; U_bot; m_fu; U_fu; m_pin]   (pressurant in ullage / surface / bulk, bottles,
%      fuel-tank ullage gas, total pressurant into the LOX tank)
x0 = [M0; U0; 1e-4; T0l; 1e-4; T0l; p.T_ln1_0; p.T_amb; p.T_amb; 0; 0; Md0; Ud0; 0; Tw0; ...
      0; 0; 0; rb*p.V_bot; rb*p.V_bot*ub; rf*p.V_fu; rf*p.V_fu*uf; 0];
iP = 14 + p.N_wall;                            % pressurant states start after the wall bands

%% 4. PROCEDURE ------------------------------------------------------------
p.t_wall0 = tic;                               % run timer for the progress lines
fprintf('\nSolving (a progress line every %g s of real time; set p.progress_dt = 0 to silence)\n', p.progress_dt);
absTol = [1e-6; 10; 1e-5; 1e-3; 1e-5; 1e-3; 1e-3; 1e-3; 1e-3; 1e-5; 1e-5; 1e-3; 10; 1e-5; repmat(1e-3,p.N_wall,1); ...
          1e-6; 1e-6; 1e-6; 1e-5; 100; 1e-6; 10; 1e-6];
opts  = odeset('RelTol',1e-5,'AbsTol',absTol,'MaxStep',5);
optsA = odeset(opts,'Events',@(t,x) ev_fill_full(t,x,p));
optsC = odeset(opts,'Events',@(t,x) ev_liquid_gone(t,x,p));

t = 0;  x = x0;  T_all = t;  X_all = x.';  ph_all = double(p.t_prechill <= 0);
plog = cell(0,10);
if p.t_prechill > 0
    [tt,xx] = ode15s(@(tt,xx) rhs2(tt,xx,p,false,true,false), [0 p.t_prechill], x, ...
        prog(opts,'pre-chill',p,false,true,false));
    [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tt,xx,0);  t = tt(end);  x = xx(end,:).';
    plog(end+1,:) = phase_line('pre-chill', t, x, p);
    x(7) = p.T_ln1_0;                         % tube 1 full of the new dewar's liquid
end
tv = p.t_vent(:);
for k = 1:p.n_fill
    t_fill0 = t;
    [tt,xx,~,~,ie] = ode15s(@(tt,xx) rhs2(tt,xx,p,true,false,false), [t t+p.tA_max], x, ...
        prog(optsA,sprintf('fill %d',k),p,true,false,false));
    [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tt,xx,1);  t = tt(end);  x = xx(end,:).';
    plog(end+1,:) = phase_line(sprintf('fill %d', k), t, x, p); %#ok<AGROW>
    if isempty(ie)
        warning('boiloff_v3:fill','Fill %d did not reach the dip tube within tA_max.', k);  break
    end
    te = p.tA_extra(min(k,numel(p.tA_extra))) + max(p.tA_min(min(k,numel(p.tA_min))) - (t - t_fill0), 0);
    if te > 0
        [tt,xx] = ode15s(@(tt,xx) rhs2(tt,xx,p,true,false,false), [t t+te], x, ...
            prog(opts,sprintf('fill %d +',k),p,true,false,false));
        [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tt,xx,1);  t = tt(end);  x = xx(end,:).';
        plog(end+1,:) = phase_line(sprintf('fill %d +', k), t, x, p); %#ok<AGROW>
    end
    if k < p.n_fill
        [tt,xx] = ode15s(@(tt,xx) rhs2(tt,xx,p,false,false,false), [t t+tv(min(k,numel(tv)))], x, ...
            prog(opts,sprintf('vent %d',k),p,false,false,false));
        [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tt,xx,2);  t = tt(end);  x = xx(end,:).';
        plog(end+1,:) = phase_line(sprintf('vent %d', k), t, x, p); %#ok<AGROW>
    elseif p.final_vent && ev_liquid_gone(t, x, p) > 0
        [tt,xx] = ode15s(@(tt,xx) rhs2(tt,xx,p,false,false,false), [t t+p.tC_max], x, ...
            prog(optsC,'final vent',p,false,false,false));
        [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tt,xx,3);  t = tt(end);  x = xx(end,:).';
        plog(end+1,:) = phase_line('final vent', t, x, p); %#ok<AGROW>
    end
end
pH = p;  pH.A_v = 0;                          % vent closed
if p.t_hold > 0
    [tt,xx] = ode15s(@(tt,xx) rhs2(tt,xx,pH,false,false,false), [t t+p.t_hold], x, ...
        prog(opts,'sealed hold',pH,false,false,false));
    [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tt,xx,4);  t = tt(end);  x = xx(end,:).';
    plog(end+1,:) = phase_line('sealed', t, x, pH);
end
t_p0 = t;  k_p0 = numel(T_all);                % pressurization starts here
if p.t_press > 0                               % solenoid open, vent and fill closed, until fire
    [tt,xx] = ode15s(@(tt,xx) rhs2(tt,xx,pH,false,false,true), [t t+p.t_press], x, ...
        prog(opts,'pressurize',pH,false,false,true));
    [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tt,xx,5);  t = tt(end);  x = xx(end,:).';
    plog(end+1,:) = phase_line('pressurized', t, x, pH);
end

%% 5. POST-PROCESSING ------------------------------------------------------
n = numel(T_all);
f = {'P','Tv','Ti','Ts','Tb','Vl','Vsw','ml','mdot_in','mdot_vent','mdot_vent_liq','m_ev','Qwl','Qwg','P_dewar', ...
     'f_vent','gam_vent','Pc','Pp','Tsw','P_fuel','P_bot','P_reg','m_sol','m_pox','m_pfu','m_abs','T_reg', ...
     'm_liq_cryo','Vv'};
for j = 1:numel(f), D.(f{j}) = zeros(n,1); end
fprintf('Post-processing %d stored points...\n', n);  t_pp = tic;  t_last = 0;
for i = 1:n
    if p.progress_dt > 0 && toc(t_pp) - t_last >= p.progress_dt
        t_last = toc(t_pp);  fprintf('  post-processing %5.1f %%  (%d / %d points, %s)\n', 100*i/n, i, n, fmt_time(t_last));
    end
    if ph_all(i) >= 4, pp = pH; else, pp = p; end
    [~,d] = rhs2(T_all(i), X_all(i,:).', pp, ph_all(i)==1, ph_all(i)==0, ph_all(i)==5);
    for j = 1:numel(f), D.(f{j})(i) = d.(f{j}); end
end
Mtank = X_all(:,1) + X_all(:,3) + X_all(:,5) + X_all(:,iP+1);     % everything in the tank, incl. pressurant
Mpress = X_all(:,iP+1) + X_all(:,iP+2) + X_all(:,iP+3);          % pressurant in the tank (ullage + dissolved)
if p.same_gas, Mcryo = Mtank; else, Mcryo = Mtank - Mpress; end    % cryogen only
if min(X_all(:,3)) < -1e-3 || min(X_all(:,5)) < -1e-3
    warning('boiloff_v3:mass','A liquid node went negative (surface %.3f, bulk %.3f kg).', min(X_all(:,3)), min(X_all(:,5)));
end

fprintf('\n%-11s %7s %8s %7s %7s %7s %8s %8s %8s %8s\n','phase end','t [s]','P [psig]','Tv [K]','Ts [K]', ...
    'Tb [K]','M [kg]','liq [kg]','in [kg]','vent [kg]');
for i = 1:size(plog,1), fprintf('%-11s %7.0f %8.1f %7.1f %7.1f %7.1f %8.2f %8.2f %8.2f %8.2f\n', plog{i,:}); end
fprintf('\n================ RESULT (stratified) ================\n');
fprintf('Final pressure                 : %.1f psig\n', (D.P(end)-p.P_amb)/psi2Pa);
fprintf('Final T: vapor / surface / bulk: %.1f / %.1f / %.1f K  (Tsat at the tank pressure %.1f K, interface %.1f K)\n', ...
    D.Tv(end), D.Ts(end), D.Tb(end), D.Tsw(end), D.Ti(end));
fprintf('%s in tank at end              : %.2f kg  (liquid %.2f kg)\n', p.cryo.label, Mcryo(end), D.m_liq_cryo(end));
if ~isnan(p.meas_mass)
    fprintf('  measured: %.2f kg (expelled in the burn; the tank held that plus residual)\n', p.meas_mass);
end
fprintf('%s loaded (net, incl. pre-chill): %.2f kg\n', p.cryo.label, X_all(end,11));
fprintf('%s vented                      : %.2f kg\n', p.cryo.label, X_all(end,10));
fprintf('Peak tank pressure             : %.1f psig at t = %.0f s\n', (max(D.P)-p.P_amb)/psi2Pa, T_all(find(D.P==max(D.P),1)));
iv = D.mdot_vent > 1e-4 & ph_all < 4;
if any(iv)
    fprintf('Vent gas exponent (venting)    : %.3f - %.3f\n', min(D.gam_vent(iv)), max(D.gam_vent(iv)));
    fprintf('Vent friction factor (venting) : %.4f - %.4f\n', min(D.f_vent(iv)), max(D.f_vent(iv)));
end
fprintf('Mass balance residual          : %.2e kg (should be ~0)\n', Mtank(end) - M0 - 2e-4 - X_all(end,11) + X_all(end,10) - X_all(end,iP+8));
if p.t_press > 0
    k0 = k_p0;  Pa = p.P_amb;
    fprintf('\n---- Pressurization: %.0f s with %s ----\n', p.t_press, p.press_gas);
    fprintf('%s tank pressure              : %.1f -> %.1f psig  (fuel tank %.1f -> %.1f psig)\n', p.cryo.label, ...
        (D.P(k0)-Pa)/psi2Pa, (D.P(end)-Pa)/psi2Pa, (D.P_fuel(k0)-Pa)/psi2Pa, (D.P_fuel(end)-Pa)/psi2Pa);
    fprintf('Bottle pressure                : %.0f -> %.0f psig  (%.2f kg %s used: %.2f kg to %s tank, %.2f kg to fuel tank)\n', ...
        (D.P_bot(k0)-Pa)/psi2Pa, (D.P_bot(end)-Pa)/psi2Pa, X_all(k0,iP+4) - X_all(end,iP+4), p.press_gas, ...
        X_all(end,iP+8) - X_all(k0,iP+8), p.cryo.label, X_all(end,iP+6) - X_all(k0,iP+6));
    fprintf('%s LIQUID (what you fire)      : %.3f -> %.3f kg  (%+.3f kg)\n', p.cryo.label, ...
        D.m_liq_cryo(k0), D.m_liq_cryo(end), D.m_liq_cryo(end) - D.m_liq_cryo(k0));
    fprintf('%s vapor in the ullage         : %.3f -> %.3f kg\n', p.cryo.label, X_all(k0,1), X_all(end,1));
    if ~p.same_gas
        fprintf('%s in the %s tank              : %.3f kg in the ullage + %.3f kg dissolved in the liquid\n', ...
            p.press_gas, p.cryo.label, X_all(end,iP+1), X_all(end,iP+2) + X_all(end,iP+3));
        fprintf('Liquid composition at fire     : %.2f %% %s by mass (surface layer %.1f %%)\n', ...
            100*(X_all(end,iP+2) + X_all(end,iP+3))/max(D.ml(end),1e-9), p.press_gas, ...
            100*X_all(end,iP+2)/max(X_all(end,3),1e-9));
    end
    fprintf('Liquid T (surface / bulk)      : %.1f / %.1f K -> %.1f / %.1f K\n', D.Ts(k0), D.Tb(k0), D.Ts(end), D.Tb(end));
    fprintf('Tsat before pressurization     : %.1f K at %.1f psig  (surface %.1f K, bulk %.1f K below it)\n', ...
        D.Tsw(k0), (D.P(k0)-Pa)/psi2Pa, D.Tsw(k0) - D.Ts(k0), D.Tsw(k0) - D.Tb(k0));
    if D.P(end) < p.tab.P(end)
        fprintf('Bulk subcooling at fire        : %.1f K below Tsat(P) = %.1f K\n', D.Tsw(end) - D.Tb(end), D.Tsw(end));
    else
        fprintf('Tank pressure is above %.0f psig (Tsat table top); %s critical pressure is %.0f psig:\n', ...
            (p.tab.P(end)-Pa)/psi2Pa, p.cryo.label, (p.P_crit-Pa)/psi2Pa);
        fprintf('  no boiling, the liquid is compressed liquid (bulk %.1f K)\n', D.Tb(end));
    end
    fprintf('Liquid volume                  : %.2f -> %.2f L  (ullage %.2f -> %.2f L)\n', 1e3*D.Vl(k0), 1e3*D.Vl(end), ...
        1e3*D.Vv(k0), 1e3*D.Vv(end));
end
fprintf('Total run time                 : %s\n', fmt_time(toc(p.t_wall0)));
fprintf('=====================================================\n');

results = struct('t',T_all,'x',X_all,'phase',ph_all,'diag',D,'M_tank',Mtank,'M_cryo',Mcryo,'M_press',Mpress, ...
                 't_press0',t_p0,'log',{plog},'p',p);
save(fullfile(here, sprintf('%s_chill_results_v3.mat', lower(p.cryo.name))),'results');

if p.make_plots
    tm = T_all/60;  iw = round([1, p.N_wall/2, p.N_wall]);
    figure('Name',sprintf('%s chill-in, %s tank, stratified (v3)', p.cryo.name, p.mw.name),'Position',[60 60 1350 800]);
    subplot(2,3,1);
    plot(tm, D.Tv, tm, D.Ts, tm, D.Tb, 'LineWidth',1.4); hold on; plot(tm, D.Tsw, 'w:', 'LineWidth',1);
    plot(tm, X_all(:,14+iw), 'LineWidth',0.9);
    xlabel('time [min]'); ylabel('T [K]'); grid on; ylim([floor(Tsat_of_P(p.P_amb, p.tab)) - 5, max(p.T_amb, 300)]);
    legend({'vapor','surface layer','bulk','T_{sat}(P)','wall bottom','wall middle','wall top'},'Location','best');
    subplot(2,3,2);
    plot(tm, (D.P-p.P_amb)/psi2Pa, 'LineWidth',1.5, 'DisplayName','tank (v3)'); hold on;
    plot(tm, (D.P_dewar-p.P_amb)/psi2Pa, '--w', 'DisplayName','dewar');
    mf = p.meas_file;  if ~isfile(mf), mf = fullfile(here, p.meas_file); end
    if p.show_meas && isfile(mf)
        Mm = readmatrix(mf);  tmeas = Mm(:,1) - (p.meas_t_fill - p.t_prechill);
        km = tmeas >= 0 & tmeas <= T_all(end);
        plot(tmeas(km)/60, p.meas_cal*(Mm(km,2)-p.meas_zero), 'r--', 'LineWidth',1.0, 'DisplayName','measured');
    end
    xlabel('time [min]'); ylabel('P [psig]'); grid on; legend('Location','best');
    subplot(2,3,3);
    plot(tm, D.Vl/p.V_dip, 'LineWidth',1.5); hold on; if p.swell, plot(tm, D.Vsw/p.V_dip, 'LineWidth',1.1); end
    yline(1,'--r','dip tube');  xlabel('time [min]'); ylabel('V / V_{dip}'); grid on;
    subplot(2,3,4);
    plot(tm, Mcryo, tm, D.m_liq_cryo, tm, X_all(:,10), tm, X_all(:,11), 'LineWidth',1.4);
    xlabel('time [min]'); ylabel('mass [kg]'); grid on; legend('in tank','liquid','vented','loaded','Location','best');
    subplot(2,3,5);
    plot(tm, D.mdot_in, tm, D.mdot_vent, tm, D.mdot_vent_liq, tm, D.m_ev, 'LineWidth',1.1);
    xlabel('time [min]'); ylabel('mass flow [kg/s]'); grid on;
    legend('in','vent total','vent liquid','evap (+) / cond (-)','Location','best');
    subplot(2,3,6);
    plot(tm, D.Qwl/1e3, tm, D.Qwg/1e3, 'LineWidth',1.1);
    xlabel('time [min]'); ylabel('wall heat [kW]'); grid on; legend('wall -> liquid','wall -> vapor','Location','best');
    if p.t_press > 0, plot_press(T_all, X_all, D, k_p0, t_p0, iP, p, here); end
end
end

%% ========================= LOCAL FUNCTIONS ================================

function c = phase_line(name, t, x, p)
    [~,d] = rhs2(t, x, p, false, false, false);
    c = {name, t, (d.P - p.P_amb)/6894.757, d.Tv, d.Ts, d.Tb, x(1)+max(x(3),0)+max(x(5),0), d.ml, x(11), x(10)};
end

function v = Vswell2(x, p)
    Vl = max(x(3),0)/sat_props(x(4), p.tab).rho_l + max(x(5),0)/sat_props(x(6), p.tab).rho_l;
    if p.swell, v = Vl/(1 - min(max(x(14),0),0.9)); else, v = Vl; end
end

function [value, isterminal, direction] = ev_fill_full(~, x, p)
    value = Vswell2(x, p) - (p.V_dip - 0.5*p.A_c*p.dh_dip);  isterminal = 1;  direction = 1;
end

function [value, isterminal, direction] = ev_liquid_gone(~, x, p)
    value = Vswell2(x, p) - (p.V_dip - p.A_c*p.dh_dip);  isterminal = 1;  direction = -1;
end

function [dx, d] = rhs2(~, x, p, fill, prechill, press)
% Time derivatives of the stratified model.
% x = [m_v; U_v; m_s; T_s; m_b; T_b; T_ln1; T_ln2; T_ln3; m_vent; m_in; M_dew; U_dew; alpha; T_wall(1:N);
%      m_pv; m_ps; m_pb; m_bot; U_bot; m_fu; U_fu; m_pin]
% m_v is cryogen vapor only; m_s and m_b are liquid masses INCLUDING dissolved pressurant (m_ps, m_pb).
% press = true while the pressurant solenoid is open.
    N  = p.N_wall;
    mv = x(1); Uv = x(2); ms = x(3); Ts = x(4); mb = x(5); Tb = x(6); T1 = x(7); T2 = x(8); T3 = x(9);
    Md = x(12); Ud = x(13); al_s = x(14); Tw = x(15:14+N);
    xp = x(15+N:end);  mpv = max(xp(1),0);  mps = max(xp(2),0);  mpb = max(xp(3),0);
    mbot = xp(4);  Ubot = xp(5);  mfu = xp(6);  Ufu = xp(7);
    tab = p.tab;  g = p.g;
    sp  = @(z) 0.5*(z + sqrt(z.^2 + 1e-10));
    sg  = @(z) 1./(1 + exp(-min(max(z,-50),50)));
    ms_e = max(ms,0);  mb_e = max(mb,0);
    ws = min(mps/max(ms_e,1e-9), 0.99);  wb = min(mpb/max(mb_e,1e-9), 0.99);   % dissolved pressurant mass fractions
    ss = sat_props(Ts, tab);  sb = sat_props(Tb, tab);
    [~, rps] = press_sat(Ts, p);  [~, rpb] = press_sat(Tb, p);                 % dissolved pressurant liquid density
    % Liquid volume: cryogen compressed to the tank pressure + dissolved pressurant. P sets V_l and V_l
    % sets P through the ullage volume, so iterate (converges in 2-3 passes; the effect is ~1%).
    P = ss.Psat;
    for it = 1:3
        ros = ss.rho_l*(1 + ss.kap*max(P - ss.Psat, 0));  rob = sb.rho_l*(1 + sb.kap*max(P - sb.Psat, 0));
        Vl  = ms_e*(1-ws)/ros + ms_e*ws/rps + mb_e*(1-wb)/rob + mb_e*wb/rpb;
        Vv  = max(p.V_tank - Vl, 1e-4);
        [Tv, Pc, Pp, uc, up] = mix_state(max(mv,1e-6), mpv, Vv, Uv, p);
        P = Pc + Pp;
    end
    rv  = max(mv,1e-6)/Vv;  hv = uc + Pc/rv;                  % cryogen vapor density, enthalpy
    rpv = mpv/Vv;                                              % pressurant density in the ullage
    if mpv > 1e-6, hpg = up + Pp/rpv; else, hpg = press_gas(Tv, P, p); hpg = hpg(6); end
    wpv = mpv/(mpv + max(mv,1e-6));                            % pressurant mass fraction of the ullage gas
    % Interface: the cryogen evaporates / condenses at the temperature where its vapor pressure over the
    % surface liquid (Raoult: x_c*Psat) equals its partial pressure. Without pressurant this is Tsat(P).
    xs  = press_molefrac(ws, p);
    Ti  = Tsat_of_P(Pc/max(1 - xs, 0.05), tab);  si = sat_props(Ti, tab);  hfg = si.h_g - si.h_l;
    % Bubbles at the wall or in the liquid must hold the TOTAL pressure: boiling, flashing and bubble
    % survival are referenced to Tsat(P). Above ~4.2 MPa (table top, 150 K for O2) there is no boiling.
    Tsw = Tsat_of_P(P, tab);  sw = sat_props(Tsw, tab);  hfgw = sw.h_g - sw.h_l;
    ml = ms_e + mb_e;  fl = ml/(ml + 0.05);
    if p.swell, al = min(max(al_s,0),0.9); else, al = 0; end
    hl = min(Vl/p.A_c, p.L_t);

    % ---- wall bands ----
    fw = min(max((hl - p.z0)/p.dz, 0), 1);
    Lg = max(p.L_t - hl, 1e-2);  Twdry = sum(Tw.*(1-fw))/max(sum(1-fw), 1e-9);
    [kg,mug,cpg,~,bg] = gas_film_mix(0.5*(Twdry + Tv), P, wpv, p);  rhog = (max(mv,1e-6) + mpv)/Vv;
    Qw_l = zeros(N,1);  Qw_g = zeros(N,1);  Qex_w = zeros(N,1);
    for i = 1:N
        if fw(i) > 0 && ml > 1e-3
            Qw_l(i) = fw(i)*p.Ab*wall_liq_flux(Tw(i), Tb, Tsw, sb, sw, max(hl,1e-2), p)*fl;
        end
        if fw(i) < 1
            Qw_g(i) = (1-fw(i))*p.Ab*nat_conv(kg,mug,cpg,rhog,bg,Tw(i)-Tv,Lg,g)*(Tw(i) - Tv);
        end
        Qex_w(i) = p.Ab*ext_heat(Tw(i), p, p.L_t, 'v', p.eps);
    end
    Qax = zeros(N,1);
    dth = p.C_ax*diff(p.theta(Tw));                           % conduction between neighbouring bands [W]
    Qax(1:end-1) = Qax(1:end-1) + dth;  Qax(2:end) = Qax(2:end) - dth;
    Qwl = sum(Qw_l);  Qwg = sum(Qw_g);
    % Wall heat goes to whichever liquid is present (split by mass). In each node it makes vapor that
    % reaches the ullage once that liquid is near saturation at the TOTAL pressure, otherwise it heats it.
    fb   = mb_e/(ml + 1e-6);                                  % share of the liquid in the bulk
    fv_b = sg((Tb - (Tsw - p.dT_sub_vap))/0.5);  fv_s = sg((Ts - (Tsw - p.dT_sub_vap))/0.5);
    f_vap = fb*fv_b + (1-fb)*fv_s;
    m_wv_b = sp(Qwl)*fb*fv_b/hfgw;  m_wv_s = sp(Qwl)*(1-fb)*fv_s/hfgw;  m_wv = m_wv_b + m_wv_s;
    Qwl_b  = Qwl*fb*(1-fv_b);       Qwl_s  = Qwl*(1-fb)*(1-fv_s);

    % ---- fill line ----
    Qex_1 = pi*p.OD1*p.L1*ext_heat(T1, p, p.OD1, 'h', p.eps_line);
    Qex_2 = pi*p.OD2*p.L2*ext_heat(T2, p, p.OD2, 'h', p.eps_line);
    Qex_3 = pi*p.OD3*p.L3*ext_heat(T3, p, p.OD3, 'h', p.eps_line);
    if p.dewar_sag
        stD = fluid_state(Md, Ud, p.pdew);  Pd = stD.P;  sd = sat_props(stD.T, tab);
    else
        Pd = p.P_d;  sd = p.s_dl;
    end
    Q1 = pi*p.ID1*p.L1*boil_flux(T1 - sd.T, sd, p.ID1, p.C_sf, g, p.dTmin_ovr);
    wet = (ml > 0.05) || fill;
    Q2 = seg_heat2(T2, p.ID2, p.L2, wet, Tb, Tsw, sb, sw, Tv, P, wpv, p);
    Q3 = seg_heat2(T3, p.ID3, p.L3, wet, Tb, Tsw, sb, sw, Tv, P, wpv, p);
    mdot_in = 0;  m_back = 0;  md0 = 0;  H_in = 0;  Q_line_tank = 0;
    if fill
        dP0 = 2000;  dP = Pd - P;  a = min(abs(dP), dP0)/dP0;
        if dP >= 0
            md0 = p.A_in*table_G(Pd - max(dP,dP0), Pd, p.G2P0, p.G2P, p.G2)*a;
        elseif p.backflow && ml > 0.05
            m_back = p.A_in*table_G(P - max(-dP,dP0), P, p.GrP0, p.GrP, p.Gr)*a;
        end
        sm  = sat_props(Tsat_of_P(0.5*(Pd + P), tab), tab);
        vl_ = 1/sm.rho_l;  vlg = 1/sm.rho_g - vl_;
        x_ex = min(max((sd.h_l - sw.h_l)/hfgw, 0), 1);
        bq  = (sp(Q1) + sp(Q2) + sp(Q3))*vlg/(2*(vl_ + 0.5*x_ex*vlg)*(sm.h_g - sm.h_l));
        mdot_in = 0.5*(-bq + sqrt(bq^2 + 4*md0^2));
        H_in = mdot_in*sd.h_l + Q1*mdot_in/(mdot_in + 0.005) + Q2 + Q3;
    else
        Q_line_tank = Q2 + Q3;
    end
    % inflow at tank pressure: flash vapor -> ullage; liquid rises to the surface if warmer than the bulk
    x_in = 0;  if mdot_in > 0, x_in = min(max((H_in/mdot_in - sw.h_l)/hfgw, 0), 1); end
    m_fl0 = mdot_in*x_in;  m_inl = mdot_in - m_fl0;
    m_fl  = m_fl0*f_vap;  m_flc = m_fl0 - m_fl;           % flash bubbles collapse in a subcooled bulk
    if x_in > 0 || mdot_in <= 0
        h_inl = sw.h_l;  T_inl = Tsw;
    else
        h_inl = H_in/mdot_in;  T_inl = interp1(tab.h_l_col, tab.T, min(max(h_inl, tab.h_l_col(1)), tab.h_l_col(end)));
    end
    w_up = sg((T_inl - Tb)/0.3)*mb_e/(mb_e + p.m_s0);       % empty tank: liquid collects at the bottom

    % ---- interface (cryogen evaporation / condensation) ----
    A_i  = p.A_c*fl;
    h_li = p.h_cond + (p.h_evap - p.h_cond)*sg((Ts - Ti)/0.2);
    Q_li = h_li*A_i*(Ts - Ti);  if Ts > Ti, Q_li = Q_li*ms_e/(ms_e + 0.01); end
    Q_vi = p.h_vi*A_i*(Tv - Ti);
    m_ev = (Q_li + Q_vi)/hfg;                                  % + evaporation, - condensation

    % ---- pressurant absorption into the surface layer ----
    % Equilibrium (Raoult) mole fraction x_eq = P_press/Psat_press(T_s), capped at 1; the rate is set by a
    % liquid-side mass-transfer coefficient beta_abs. Negative = the dissolved gas comes back out.
    m_abs = 0;  Q_abs = 0;
    if p.press_dissolves && ml > 0.05
        [Psp, ~, hpl] = press_sat(Ts, p);
        w_eq  = press_massfrac(min(Pp/Psp, 1), p);
        m_abs = p.beta_abs*A_i*ss.rho_l*(w_eq - ws);
        if m_abs > 0, m_abs = m_abs*mpv/(mpv + 1e-4); else, m_abs = m_abs*mps/(mps + 1e-4); end
        Q_abs = m_abs*(hpg - hpl);                            % condensation heat released in the surface layer
    end

    % ---- surface layer <-> bulk ----
    h_sb = p.h_sb_stable + (p.h_sb_mix - p.h_sb_stable)*sg((Tb - Ts)/0.3);
    Q_sb = h_sb*p.A_c*fl*(Ts - Tb);
    if ms_e < p.m_s0, m_x = (p.m_s0 - ms_e)/p.tau_up*mb_e/(mb_e + 0.05);   % bulk -> surface
    else,             m_x = -(ms_e - p.m_s0)/p.tau_down;                  % surface -> bulk
    end
    % dissolved pressurant into the surface: carried by m_x, plus slow exchange (N2-rich liquid is lighter)
    if m_x > 0, m_px = m_x*wb; else, m_px = m_x*ws; end
    m_px = m_px + p.beta_sb*p.A_c*fl*sb.rho_l*(wb - ws);
    m_boil = sp(Tb - Tsw)*mb_e*sb.cp_l/(hfgw*p.tau_boil);

    % ---- vent (ullage gas; the pressurant leaves in proportion to its mass fraction) ----
    rho_vent = (max(mv,1e-6) + mpv)/Vv;
    [m_vap, f_vt, g_vt] = vent_flow(P, rho_vent, Tv, p);  m_liq = 0;  w = 0;  Vsw = Vl;
    if ml > 0.05
        Vsw = Vl/(1 - al);  w = sg((Vsw/p.A_c - p.h_dip)/(p.dh_dip/4));
        rho_m = al*sw.rho_g + (1-al)*ss.rho_l;
        m_liq = p.A_v*max(interp1(p.vl_P, p.vl_G, P, 'linear', 'extrap'), 0)*sqrt(rho_m/ss.rho_l);
    end
    m_vl = w*m_liq*ms_e/(ms_e + 0.01);
    m_vv = (1-w)*m_vap;  m_vv_c = m_vv*(1 - wpv);  m_vv_p = m_vv*wpv;
    m_src = 0;  if prechill, m_src = vent_flow(p.P_hold, rv*p.P_hold/P, Tv, p); end

    % ---- pressurant supply: bottles -> regulator -> solenoid -> manifold -> check valves -> tanks ----
    [Tbot, Pbot, ubot] = press_state(mbot, p.V_bot, Ubot, p);
    [Tfu,  Pfu]        = press_state(mfu,  p.V_fu,  Ufu,  p);
    h0 = ubot + Pbot*p.V_bot/mbot;                            % bottle gas enthalpy (kept through the throttles)
    m_sol = 0;  m_pox = 0;  m_pfu = 0;  Qbot = 0;  Qfu = 0;  Pm = NaN;  Preg = NaN;  T0 = NaN;
    if press
        T0 = press_T_of_hP(h0, min(p.P_reg, Pbot), p);
        [Pm, m_sol, m_pox, m_pfu, Preg] = press_network(Pbot, T0, P, Pfu, p);
        gb = press_gas(0.5*(Tbot + p.T_amb), Pbot, p);
        Qbot = p.A_bot*nat_conv(gb(1),gb(2),gb(3),gb(4),gb(5),p.T_amb - Tbot,p.L_bot,g)*(p.T_amb - Tbot);
        gf = press_gas(0.5*(Tfu + p.T_amb), Pfu, p);  Lfu = p.V_fu/p.A_fuc;
        Qfu = (pi*p.D_fu*Lfu + 2*p.A_fuc)*nat_conv(gf(1),gf(2),gf(3),gf(4),gf(5),p.T_amb - Tfu,Lfu,g)*(p.T_amb - Tfu);
    end
    if p.same_gas, m_pc = m_pox; m_pp = 0; else, m_pc = 0; m_pp = m_pox; end   % LN2 + GN2: one species

    % ---- mass ----
    dmv = m_fl + m_ev + m_wv + m_boil + m_src - m_vv_c + m_pc;
    dmpv = m_pp - m_abs - m_vv_p;
    m_back_b = m_back*fb;  m_back_s = m_back - m_back_b;      % backflow / collapsed bubbles split by mass
    m_flc_b  = m_flc*fb;   m_flc_s  = m_flc - m_flc_b;
    dms = m_inl*w_up + m_flc_s - m_ev + m_x - m_vl - m_wv_s - m_back_s + m_abs;
    dmb = m_inl*(1-w_up) + m_flc_b - m_x - m_wv_b - m_boil - m_back_b;
    dmps = m_abs + m_px - (m_vl + m_back_s)*ws;
    dmpb = -m_px - m_back_b*wb;
    % ---- liquid temperatures (inflows at h_i change T; outflows at h(T) do not) ----
    cs = max(ms_e,0.05)*ss.cp_l;  cb = max(mb_e,0.05)*sb.cp_l;
    dTs = m_inl*w_up*(h_inl - ss.h_l) + m_flc_s*(sw.h_g - ss.h_l) - m_wv_s*(sw.h_l - ss.h_l) + Qwl_s ...
          - m_ev*(si.h_l - ss.h_l) - Q_li - Q_sb + Q_abs;
    if m_x > 0, dTs = dTs + m_x*(sb.h_l - ss.h_l); end
    dTs = dTs/cs;
    dTb = m_inl*(1-w_up)*(h_inl - sb.h_l) + m_flc_b*(sw.h_g - sb.h_l) - m_wv_b*(sw.h_l - sb.h_l) - m_boil*(sw.h_g - sb.h_l) ...
          + Qwl_b + Q_sb + Q_line_tank*fl;
    if m_x < 0, dTb = dTb - m_x*(ss.h_l - sb.h_l); end
    dTb = dTb/cb;
    % ---- ullage energy (P dV from the change in liquid volume; the ~1% compression of the liquid is
    %      included in V_l but its tiny P dV work is not) ----
    drs = (sat_props(Ts+0.05,tab).rho_l - sat_props(Ts-0.05,tab).rho_l)/0.1;
    drb = (sat_props(Tb+0.05,tab).rho_l - sat_props(Tb-0.05,tab).rho_l)/0.1;
    dVl = (dms - dmps)/ros + dmps/rps - ms_e*(1-ws)/ros^2*drs*dTs ...
        + (dmb - dmpb)/rob + dmpb/rpb - mb_e*(1-wb)/rob^2*drb*dTb;
    dUv = (m_fl + m_ev + m_wv + m_boil)*sw.h_g + m_src*p.h_src - m_vv_c*hv - m_vv_p*hpg ...
          + Qwg - Q_vi + Q_line_tank*(1 - fl) + P*dVl + m_pox*h0 - m_abs*hpg;

    % ---- wall, line, void ----
    dTwall = (Qex_w - Qw_l - Qw_g + Qax)./(p.mb*cp_mat(Tw, p.mw));
    dT1 = (Qex_1 - Q1)/(p.m_ln1*cp_mat(T1, p.mln));
    dT2 = (Qex_2 - Q2)/(p.m_ln2*cp_mat(T2, p.mln));
    if p.L3 > 0, dT3 = (Qex_3 - Q3)/(p.m_ln3*cp_mat(T3, p.mln)); else, dT3 = 0; end
    dal = -al_s/5;
    if p.swell && ml > 0.05
        j   = (0.5*m_wv + 0.5*m_boil + m_fl)/(sw.rho_g*p.A_c);
        Ugj = p.Ugj_coef*(sw.sigma*g*(sw.rho_l - sw.rho_g)/sw.rho_l^2)^0.25;
        dal = (j/(p.C0*j + Ugj) - al_s)/(max(hl,0.02)/(p.C0*j + Ugj));
    end
    dx = [dmv; dUv; dms; dTs; dmb; dTb; dT1; dT2; dT3; m_vv + m_vl; mdot_in - m_back + m_src; ...
          -mdot_in + m_back; -mdot_in*sd.h_l + m_back_b*sb.h_l + m_back_s*ss.h_l; dal; dTwall; ...
          dmpv; dmps; dmpb; -m_sol; -m_sol*h0 + Qbot; m_pfu; m_pfu*h0 + Qfu; m_pox];
    d = struct('P',P,'Tv',Tv,'Ti',Ti,'Ts',Ts,'Tb',Tb,'Vl',Vl,'Vsw',Vsw,'ml',ml, ...
               'mdot_in',mdot_in - m_back + m_src,'mdot_vent',m_vv + m_vl,'mdot_vent_liq',m_vl, ...
               'm_ev',m_ev,'Qwl',Qwl,'Qwg',Qwg,'P_dewar',Pd,'f_vent',f_vt,'gam_vent',g_vt, ...
               'Pc',Pc,'Pp',Pp,'Tsw',Tsw,'P_fuel',Pfu,'P_bot',Pbot,'P_reg',Preg,'P_man',Pm, ...
               'm_sol',m_sol,'m_pox',m_pox,'m_pfu',m_pfu,'m_abs',m_abs,'T_reg',T0,'T_bot',Tbot,'T_fuel',Tfu, ...
               'm_liq_cryo',ml - mps - mpb,'Vv',Vv);
end

function q = wall_liq_flux(Tw, Tl, Tsw, sl, sw, L, p)
% Wall -> liquid heat flux [W/m^2]. Single-phase natural convection to the liquid at Tl, or boiling
% (full boiling curve vs T_wall - Tsat at the total pressure), whichever is larger. When the tank is
% pressurized the liquid is subcooled and the wall is below Tsat(P): no boiling, convection only.
    qc = nat_conv(sl.k_l, sl.mu_l, sl.cp_l, sl.rho_l, sl.beta_l, Tw - Tl, L, p.g)*(Tw - Tl);
    qb = 0;
    if Tw > Tsw, qb = boil_flux(Tw - Tsw, sw, L, p.C_sf, p.g, p.dTmin_ovr); end
    q = max(qc, qb);
end

function Q = seg_heat2(Tw, ID, L, wetted, Tl, Tsw, sl, sw, Tv, P, wpv, p)
% Heat from an uninsulated fill-line segment into the fluid inside it.
    if L <= 0, Q = 0; return; end
    A = pi*ID*L;
    if wetted
        Q = A*wall_liq_flux(Tw, Tl, Tsw, sl, sw, ID, p);
    else
        [k,mu,cp,rho,b] = gas_film_mix(0.5*(Tw + Tv), P, wpv, p);
        Q = A*nat_conv(k,mu,cp,rho,b,Tw-Tv,ID,p.g)*(Tw - Tv);
    end
end

function [T, Pc, Pp, uc, up] = mix_state(mc, mp, V, U, p)
% Ullage gas = cryogen vapor (mass mc) + pressurant (mp) sharing volume V at one temperature, with total
% internal energy U. Dalton: each gas fills V at its own density; P = Pc + Pp. Real-gas (T, rho) tables;
% at fixed density u rises with T, so T is read straight off the summed energy curve.
    tab = p.tab;
    [ucc, pcc] = vcol(tab.vlogR, tab.vapU, tab.vapP, mc/V);
    if mp > 1e-9
        [upc, ppc] = vcol(tab.plogR, tab.pvU, tab.pvP, mp/V);
    else
        upc = zeros(size(ucc));  ppc = upc;
    end
    Uc = mc*ucc + mp*upc;
    Uq = min(max(U, Uc(1)), Uc(end));
    T  = interp1(Uc, tab.vapT, Uq);
    Pc = interp1(tab.vapT, pcc, T);  Pp = interp1(tab.vapT, ppc, T);
    uc = interp1(tab.vapT, ucc, T);  up = interp1(tab.vapT, upc, T);
end

function [T, P, u] = press_state(m, V, U, p)
% Temperature and pressure of pure pressurant gas (bottles, fuel-tank ullage).
    [uc, pc] = vcol(p.tab.plogR, p.tab.pvU, p.tab.pvP, m/V);
    u = min(max(U/m, uc(1)), uc(end));
    T = interp1(uc, p.tab.vapT, u);  P = interp1(p.tab.vapT, pc, T);
end

function [uc, pc] = vcol(lR, Ut, Pt, rho)
% Columns u(T) and P(T) of a (T, log rho) table at density rho.
    lr = min(max(log(max(rho,1e-12)), lR(1)), lR(end));
    j  = find(lR <= lr, 1, 'last');  j = min(max(j,1), numel(lR)-1);
    w  = (lr - lR(j))/(lR(j+1) - lR(j));
    uc = (1-w)*Ut(:,j) + w*Ut(:,j+1);  pc = (1-w)*Pt(:,j) + w*Pt(:,j+1);
end

function v = press_gas(T, P, p)
% Pressurant gas properties [k mu cp rho beta h] at (T, P).
    T = min(max(T, p.tab.pgT(1)), p.tab.pgT(end));  P = min(max(P, p.tab.pgP(1)), p.tab.pgP(end));
    v = p.tab.pgas([T, log(P)]);  v = v(:);
end

function T = press_T_of_hP(h, P, p)
% Temperature of pressurant gas with enthalpy h at pressure P (after the regulator: throttling keeps h).
    P  = min(max(P, p.tab.pgP(1)), p.tab.pgP(end));
    hc = squeeze(p.tab.pgas({p.tab.pgTv, log(P)}));  hc = hc(:,6);
    ok = [true; diff(hc) > 0];                     % below saturation the table is flat (saturated vapor)
    hc = hc(ok);  Tc = p.tab.pgTv(ok);
    T  = interp1(hc, Tc, min(max(h, hc(1)), hc(end)));
end

function [k, mu, cp, rho, beta] = gas_film_mix(Tf, P, wp, p)
% Ullage gas film properties: cryogen-vapor and pressurant properties blended by mass fraction wp.
    [k, mu, cp, rho, beta] = gas_film(Tf, P, p);
    if wp > 1e-6
        v  = press_gas(Tf, P, p);
        k  = (1-wp)*k + wp*v(1);  mu = (1-wp)*mu + wp*v(2);  cp = (1-wp)*cp + wp*v(3);
        rho = (1-wp)*rho + wp*v(4);  beta = (1-wp)*beta + wp*v(5);
    end
end

function [Ps, rho_l, h_l] = press_sat(T, p)
% Saturated-liquid pressurant properties at T (clamped to its table). A non-condensable pressurant
% (helium) returns Psat = Inf so nothing dissolves.
    if ~p.press_dissolves && ~p.same_gas, Ps = Inf; rho_l = 1000; h_l = 0; return, end
    T = min(max(T, p.tab.psT(1)), p.tab.psT(end));
    r = p.tab.psF(T);  Ps = r(1);  rho_l = r(2);  h_l = r(3);
end

function x = press_molefrac(w, p)
% Pressurant mole fraction from its mass fraction in the liquid.
    a = w/p.M_p;  b = (1-w)/p.M_c;  x = a/(a + b + 1e-30);
end

function w = press_massfrac(x, p)
% Pressurant mass fraction from its mole fraction in the liquid.
    a = x*p.M_p;  b = (1-x)*p.M_c;  w = a/(a + b);
end

function [Pm, m_sol, m_ox, m_fu, Preg] = press_network(Pbot, T0, Pox, Pfu, p)
% Regulator (outlet droops with flow) -> solenoid orifice -> manifold -> check valve + branch -> each tank.
% The manifold volume is negligible, so its pressure is where the solenoid flow equals the branch flows.
% The manifold pressure is solved to near machine precision: a loose tolerance adds noise to the flows
% that the stiff solver sees when it estimates its Jacobian, and that also makes it take tiny steps.
    Phi = min(p.P_reg, Pbot);  Plo = min([Pox, Pfu, Phi]) - 2e5;
    f = @(Pm) net_resid(Pm, Pbot, T0, Pox, Pfu, p);
    if f(Phi) >= 0
        Pm = Phi;
    elseif f(Plo) <= 0
        Pm = Plo;
    else
        Pm = fzero(f, [Plo Phi], optimset('TolX', 1e-9));
    end
    [~, m_sol, m_ox, m_fu, Preg] = net_resid(Pm, Pbot, T0, Pox, Pfu, p);
end

function [r, m_sol, m_ox, m_fu, Preg] = net_resid(Pm, Pbot, T0, Pox, Pfu, p)
    m_sol = 0;  Preg = min(p.P_reg, Pbot);
    for k = 1:6                                                % regulator droop: P_reg = P_set - k*mdot
        Preg  = min(p.P_reg - p.k_droop*m_sol, Pbot);
        m_sol = gas_orifice(Preg, Pm, T0, p.CdA_sol, p);
    end
    m_ox = gas_orifice(Pm, Pox + p.P_crack, T0, p.CdA_br, p);
    m_fu = gas_orifice(Pm, Pfu + p.P_crack, T0, p.CdA_br, p);
    r = m_sol - m_ox - m_fu;
end

function m = gas_orifice(Pu, Pd, T0, CdA, p)
% Ideal-gas orifice (isentropic nozzle x discharge coefficient), choked or subsonic, made smooth everywhere
% so the stiff solver never meets a corner:
%   dPe = sqrt(dP^2 + dP0^2) sets the nozzle pressure ratio (dP0 = 2000 Pa = 0.3 psi), and the factor
%   dP+/dPe makes the flow linear in dP near zero instead of the real sqrt(dP) with its infinite slope.
%   dP+ = 0.5*(dP + sqrt(dP^2 + dPs^2)) is a smooth max(dP, 0): the check valve. dPs = 50 Pa (0.007 psi),
%   so the "leak" it allows when closed is negligible.
% Above a few psi of pressure drop this is the normal orifice equation to within 0.1 %.
    dP0 = 2000;  dPs = 50;  gm = p.gam_p;  dP = Pu - Pd;
    dPp = 0.5*(dP + sqrt(dP^2 + dPs^2));
    dPe = sqrt(dP^2 + dP0^2);  r = Pd/(Pd + dPe);
    rc = (2/(gm+1))^(gm/(gm-1));
    if r <= rc
        phi = (2/(gm+1))^((gm+1)/(2*(gm-1)));
    else
        phi = sqrt(2/(gm-1)*(r^(2/gm) - r^((gm+1)/gm)));
    end
    m = CdA*phi*(Pd + dPe)*sqrt(gm/(p.R_p*T0))*dPp/dPe;
end

function plot_press(T_all, X_all, D, k0, t0, iP, p, here)
% Pressurization figure: model vs the measured traces, time from solenoid open.
    psi = 6894.757;  Pa = p.P_amb;  k = k0:numel(T_all);  tt = T_all(k) - t0;
    M = [];
    mf = p.meas_press_file;  if ~isempty(mf) && ~isfile(mf), mf = fullfile(here, mf); end
    if ~isempty(mf) && isfile(mf), M = readmatrix(mf, 'NumHeaderLines', 4); end
    figure('Name', sprintf('%s pressurization (v3)', p.cryo.name), 'Position', [90 90 1300 760]);
    subplot(2,2,1);
    plot(tt, (D.P(k)-Pa)/psi, tt, (D.P_fuel(k)-Pa)/psi, tt, (D.P_reg(k)-Pa)/psi, 'LineWidth', 1.4); hold on;
    lg = {sprintf('%s tank', p.cryo.label), 'fuel tank', 'regulator outlet'};
    if ~isempty(M)
        plot(M(:,1), M(:,2), 'r--', M(:,1), M(:,3), 'r-.', M(:,1), M(:,5), 'r:', 'LineWidth', 1.0);   % measured: red
        lg = [lg, {'measured LOX', 'measured fuel', 'measured regulator'}];
    end
    xlabel('time from solenoid open [s]'); ylabel('P [psig]'); grid on; legend(lg, 'Location', 'southeast');
    subplot(2,2,2);
    plot(tt, (D.P_bot(k)-Pa)/psi, 'LineWidth', 1.4); hold on;
    if ~isempty(M), plot(M(:,1), M(:,4), 'r--', 'LineWidth', 1.0); legend('model','measured','Location','best'); end
    xlabel('time from solenoid open [s]'); ylabel('bottle P [psig]'); grid on;
    subplot(2,2,3);
    plot(tt, D.Tv(k), tt, D.Ts(k), tt, D.Tb(k), tt, D.Tsw(k), 'w:', tt, D.Ti(k), 'w--', 'LineWidth', 1.2);
    xlabel('time from solenoid open [s]'); ylabel('T [K]'); grid on;
    legend('ullage gas','surface layer','bulk liquid','T_{sat}(P_{total})','interface (O_2 partial P)','Location','best');
    subplot(2,2,4);
    m0 = D.m_liq_cryo(k0);
    plot(tt, D.m_liq_cryo(k) - m0, tt, X_all(k,1) - X_all(k0,1), tt, X_all(k,iP+1), ...
         tt, X_all(k,iP+2) + X_all(k,iP+3), 'LineWidth', 1.3);
    xlabel('time from solenoid open [s]'); ylabel('mass [kg]'); grid on;
    legend(sprintf('%s liquid change', p.cryo.label), sprintf('%s vapor change', p.cryo.label), ...
           sprintf('%s in ullage', p.press_gas), sprintf('%s dissolved', p.press_gas), 'Location', 'best');
end

function o = prog(opts, name, p, fill, prechill, press)
% Adds the progress printer to an odeset structure for one phase.
    if p.progress_dt > 0
        o = odeset(opts, 'OutputFcn', @(t, y, flag) progress_out(t, y, flag, name, p, fill, prechill, press));
    else
        o = opts;
    end
end

function status = progress_out(t, y, flag, name, p, fill, prechill, press)
% ode15s output function: called after every accepted step. Prints a line every p.progress_dt seconds of
% real time with where the solver is in the phase, how big its steps are, and the tank state, so a slow
% run (steps of ~0.1-5 s) can be told apart from a stuck one (steps of ~1e-4 s or less, sim time frozen).
    persistent t0 t1 tw0 tlast nstep tprev tsim_last name0
    status = false;
    switch flag
        case 'init'
            t0 = t(1);  t1 = t(end);  tw0 = toc(p.t_wall0);  tlast = tw0;  nstep = 0;  tprev = t0;
            tsim_last = t0;  name0 = name;
            fprintf('[%s] %-11s start  t = %.0f s, runs to at most %.0f s\n', fmt_time(tw0), name, t0, t1);
        case ''
            nstep = nstep + numel(t);
            dt = t(end) - tprev;  if numel(t) > 1, dt = (t(end) - tprev)/numel(t); end
            tprev = t(end);
            tw = toc(p.t_wall0);
            if tw - tlast >= p.progress_dt
                [~, d] = rhs2(t(end), y(:,end), p, fill, prechill, press);
                rate = (t(end) - tsim_last)/(tw - tlast);           % sim seconds per real second
                tlast = tw;  tsim_last = t(end);
                pct = 100*(t(end) - t0)/max(t1 - t0, eps);
                msg = sprintf('[%s] %-11s t = %7.1f s (%4.1f %% of span) | %5d steps, step %.2g s, %.2g sim s per s | P %6.1f psig', ...
                    fmt_time(tw), name0, t(end), pct, nstep, dt, rate, (d.P - p.P_amb)/6894.757);
                if fill
                    msg = [msg sprintf(', fill %.3f kg/s, level %.1f %% of dip tube', d.mdot_in, 100*d.Vsw/p.V_dip)]; %#ok<AGROW>
                else
                    msg = [msg sprintf(', liquid %.2f kg', d.m_liq_cryo)]; %#ok<AGROW>
                end
                if press
                    msg = [msg sprintf(', fuel %.1f psig, bottles %.0f psig, N2 in %.1f g/s', ...
                        (d.P_fuel - p.P_amb)/6894.757, (d.P_bot - p.P_amb)/6894.757, 1e3*d.m_pox)]; %#ok<AGROW>
                end
                if dt < 1e-4, msg = [msg '  << tiny steps: solver struggling here']; end %#ok<AGROW>
                fprintf('%s\n', msg);
            end
        case 'done'
            tw = toc(p.t_wall0);  fprintf('[%s] %-11s done   %d steps in %s\n', fmt_time(tw), name0, nstep, fmt_time(tw - tw0));
    end
end

function s = fmt_time(sec)
% Elapsed time as m:ss (or h:mm:ss).
    sec = round(sec);  h = floor(sec/3600);  m = floor(mod(sec,3600)/60);  s2 = mod(sec,60);
    if h > 0, s = sprintf('%d:%02d:%02d', h, m, s2); else, s = sprintf('%d:%02d', m, s2); end
end

function [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tt,xx,code)
    T_all  = [T_all; tt(2:end)];  X_all = [X_all; xx(2:end,:)];
    ph_all = [ph_all; code*ones(numel(tt)-1,1)];
end

function st = fluid_state(M, U, p)
% FLUID_STATE  Tank (or dewar) state from total mass M and internal energy U.
%   Two-phase: solve T so that the volume and energy constraints both hold.
%   Single-phase vapor (warm start, no liquid): real-gas lookup from density and energy.
    tab = p.tab;
    M = max(M, 1e-6);  rb = M/p.V_tank;  ub = U/M;
    st = struct('two_phase',false,'T',NaN,'Tg',NaN,'P',NaN,'ml',0,'mg',M, ...
                'Vl',0,'rho_l',NaN,'s',[],'rho',rb);
    single = (rb <= tab.rho_g(1));
    if ~single
        Tdew   = tab.TofRhoG(min(rb, tab.rho_g(end)));
        single = (ub >= tab.ug(Tdew));                       % more energy than saturated vapor
    end
    if single
        try
            Tg = cprops('T','D',rb,'U',ub,p.fluid);
            st.Tg = Tg;  st.T = Tg;  st.P = cprops('P','D',rb,'T',Tg,p.fluid);
        catch err
            error('boiloff_v3:state','Single-phase state lookup failed (rho = %.3g kg/m3, u = %.4g J/kg): %s', ...
                  rb, ub, err.message);
        end
        return
    end
    % Two-phase: liquid exists below the dew point and below the bubble point.
    rbc  = min(max(rb, tab.rho_l_asc(1)), tab.rho_l_asc(end));
    Tbub = tab.TofRhoL(rbc);
    Tlo = tab.T(1);  Thi = min(Tdew, Tbub);
    fres = @(T) twophase_resid(T, M, U, p);
    flo = fres(Tlo);  fhi = fres(Thi);
    if flo*fhi > 0
        if abs(flo) < abs(fhi), T = Tlo; else, T = Thi; end
    else
        T = fzero(fres, [Tlo Thi], optimset('TolX',1e-9));
    end
    s  = sat_props(T, tab);
    mg = (p.V_tank - M/s.rho_l)/(1/s.rho_g - 1/s.rho_l);
    mg = min(max(mg,0), M);  ml = M - mg;
    st.two_phase = true;  st.T = T;  st.Tg = T;  st.P = s.Psat;
    st.ml = ml;  st.mg = mg;  st.rho_l = s.rho_l;  st.Vl = ml/s.rho_l;  st.s = s;
end

function f = twophase_resid(T, M, U, p)
% Energy residual m_l*u_l + m_g*u_g - U, with the masses from the volume constraint.
    r  = p.tab.Fres(T);                                % [rho_l rho_g u_l u_g]
    mg = (p.V_tank - M/r(1))/(1/r(2) - 1/r(1));
    f  = (M - mg)*r(3) + mg*r(4) - U;
end

function h = nat_conv(k, mu, cp, rho, beta, dT, L, g)
% Churchill-Chu natural convection, vertical surface.
    L  = max(L, 1e-3);
    nu = mu/rho;  alpha = k/(rho*cp);  Pr = nu/alpha;
    Ra = g*beta*sqrt(dT^2 + 0.01)*L^3/(nu*alpha);      % smooth |dT|
    Nu = (0.825 + 0.387*Ra^(1/6)/(1 + (0.492/Pr)^(9/16))^(8/27))^2;
    h  = Nu*k/L;
end

function q = boil_flux(dTe, s, Lc, Csf, g, dTmin_ovr)
% Wall heat flux [W/m^2] vs wall superheat dTe = T_wall - T_sat.
% Natural convection / nucleate (Rohsenow) / transition / film (Berenson), blended smoothly.
    hfg  = s.h_g - s.h_l;
    Prl  = s.mu_l*s.cp_l/s.k_l;
    hc   = nat_conv(s.k_l, s.mu_l, s.cp_l, s.rho_l, s.beta_l, dTe, Lc, g);
    qconv = hc*dTe;
    if dTe <= 0, q = qconv; return; end
    drho = s.rho_l - s.rho_g;
    capL = sqrt(g*drho/s.sigma);                                                    % 1/Laplace length
    qnb   = s.mu_l*hfg*capL*(s.cp_l*dTe/(Csf*hfg*Prl^1.7))^3;                       % Rohsenow
    qCHF  = 0.131*hfg*sqrt(s.rho_g)*(s.sigma*g*drho)^0.25;                          % Zuber CHF
    dTCHF = (Csf*hfg*Prl^1.7/s.cp_l)*(qCHF/(s.mu_l*hfg*capL))^(1/3);
    qmin  = 0.09*s.rho_g*hfg*(s.sigma*g*drho/(s.rho_l+s.rho_g)^2)^0.25;             % Zuber min. flux
    dTmin = 0.127*(s.rho_g*hfg/s.k_g)*(g*drho/(s.rho_l+s.rho_g))^(2/3) ...          % Berenson
            *(s.sigma/(g*drho))^0.5*(s.mu_g/(g*drho))^(1/3);
    if ~isnan(dTmin_ovr), dTmin = dTmin_ovr; end
    dTmin = max(dTmin, 1.05*dTCHF);
    Lb    = 1/capL;                                                                 % sqrt(sigma/(g*drho))
    hpfg  = hfg + 0.4*s.cp_g*dTe;
    hfb   = 0.425*(s.k_g^3*s.rho_g*drho*g*hpfg/(s.mu_g*dTe*Lb))^0.25;               % Berenson film
    qfb   = hfb*dTe;
    m     = log(qmin/qCHF)/log(dTmin/dTCHF);
    qa = (qconv^4 + qnb^4)^0.25;                         % below CHF
    qb = min(max(qCHF*(dTe/dTCHF)^m, qmin), qCHF);       % transition
    sig = @(x,x0,w) 1/(1 + exp(-(x-x0)/w));
    w1 = sig(dTe, dTCHF, max(0.05*dTCHF, 0.5));
    w2 = sig(dTe, dTmin, max(0.05*dTmin, 0.5));
    q  = (1-w1)*qa + w1*((1-w2)*qb + w2*qfb);
end

function q = ext_heat(Tw, p, L, kind, eps)
% Heat flux INTO a wall from ambient air: natural convection + radiation.
%   kind = 'v' vertical plate, 'h' horizontal cylinder (Churchill-Chu).
    Tf = min(max(0.5*(Tw + p.T_amb), p.air.T(1)), p.air.T(end));
    a  = p.air.F(Tf);                                    % [k mu cp rho]
    k = a(1); mu = a(2); cp = a(3); rho = a(4);
    nu = mu/rho;  alpha = k/(rho*cp);  Pr = nu/alpha;
    dT = p.T_amb - Tw;
    Ra = p.g*(1/Tf)*abs(dT)*L^3/(nu*alpha);
    if kind == 'v'
        Nu = (0.825 + 0.387*Ra^(1/6)/(1 + (0.492/Pr)^(9/16))^(8/27))^2;
    else
        Nu = (0.60  + 0.387*Ra^(1/6)/(1 + (0.559/Pr)^(9/16))^(8/27))^2;
    end
    h = Nu*k/L;
    q = p.q_ext_mult*(h*dT + eps*p.sigSB*(p.T_amb^4 - Tw^4));
end

function [mdot, fD, kap] = vent_flow(P, rho0, T0, p)
% Vapor vent flow through the dip tube to ambient (Fanno flow with friction).
% Phi(r, gam, K) = G/sqrt(gam*P*rho0) from setup_flows. gam is the ullage gas's isentropic exponent
% at (T0, P) and K = f*L/D + minor losses, with f solved together with the flow (f depends on
% Re = G*D/mu, and G depends on f). Fixed-point iteration; converges in a few steps.
    dP  = P - p.P_amb;
    dPe = 0.5*(dP + sqrt(dP^2 + 5^2));
    r   = min(max(p.P_amb/(p.P_amb + dPe), 0), 1);
    Tq  = min(max(T0, p.tab.gasT(1)), p.tab.gasT(end));
    Pq  = min(max(P,  p.tab.gasP(1)), p.tab.gasP(end));
    v   = p.tab.gasMuKap([Tq, log(Pq)]);  mu = v(1);  kap = v(2);
    if ~isnan(p.gam_v), kap = p.gam_v; end
    kq  = min(max(kap, p.vt_gam(1)), p.vt_gam(end));
    a0  = sqrt(kap*P*rho0);
    phiK = @(fD) p.vt_F(r, kq, min(max(log(fD*p.L_vent/p.ID_vent + p.K_vent_minor), p.vt_lK(1)), p.vt_lK(end)));
    if ~isnan(p.f_vent)
        fD = p.f_vent;
    else
        fD = p.f_vent0;
        for it = 1:3
            G  = max(phiK(fD), 0)*a0;
            fD = f_darcy(G*p.ID_vent/mu, 1.5e-6/p.ID_vent);
        end
    end
    phi  = max(phiK(fD), 0);
    mdot = p.A_v*phi*a0;
end

function G = table_G(Pb, P0, P0g, Pg, Gt)
% Mass flux [kg/m2/s] (on the 1/2 in tube area) from a source of saturated liquid at P0 into a
% back pressure Pb, from tables built at source pressures P0g. Each table is evaluated at the
% same pressure DROP (P0 - Pb) before blending, so G goes to zero as Pb -> P0.
    P0 = min(max(P0, P0g(1)), P0g(end));
    j  = min(max(find(P0g <= P0, 1, 'last'), 1), numel(P0g)-1);
    w  = (P0 - P0g(j))/(P0g(j+1) - P0g(j));
    dP = P0 - Pb;
    qa = min(max(P0g(j)   - dP, Pg(1)), Pg(end));
    qb = min(max(P0g(j+1) - dP, Pg(1)), Pg(end));
    G  = max((1-w)*interp1(Pg, Gt(j,:), qa) + w*interp1(Pg, Gt(j+1,:), qb), 0);
end

function p = setup_flows(p)
% Build the inflow table G(P_tank, P_dewar), the vapor vent table Phi(r) and the
% two-phase vent table G_liq(P_tank).
    tab = p.tab;
    % Fill path in series. Kin = loss at the segment inlet, on that segment's velocity.
    Kc  = 0.5*(1 - (p.ID1/p.ID_hose)^2);                  % hose -> 1/2 in tube contraction
    sg3 = (p.ID2/p.ID3)^2;                                % 1/2 -> 3/4 in area ratio
    K_ent = p.K_ent;
    K_dwh = 0;                                            % dewar valve K on the hose velocity
    if ~isnan(p.Cv_dewar), K_dwh = 894*(p.ID_hose/0.0254)^4/p.Cv_dewar^2; end
    sdw = (p.ID_dw/p.ID_hose)^2;                          % dewar tube -> hose area ratio
    if p.L_dw > 0
        K_hose = (1-sdw)^2/sdw^2 + K_dwh;                 % valve + expansion into the hose
    else
        K_hose = 0;  K_ent = K_ent + K_dwh;               % no dewar tube: valve at the hose inlet
    end
    % dewar tube -> valve -> hose -> 1/2 in tube -> ball valve -> 1/2 in tube -> 3/4 in tube
    segs = struct('L',  {p.L_dw,  p.L_hose,   p.L1,         p.L2,     p.L3}, ...
                  'D',  {p.ID_dw, p.ID_hose,  p.ID1,        p.ID2,    p.ID3}, ...
                  'eps',{1.5e-6,  p.eps_hose, p.eps_ln,     p.eps_ln, p.eps_ln}, ...
                  'Kloc',{0,      0,          0,            0,        0}, ...
                  'Kin',{0,       K_hose,     Kc + p.K_fit, p.K_valve, (1-sg3)^2/sg3^2});   % Borda-Carnot expansion
    segs = segs([segs.L] > 0);
    fr = [0.5 0.6 0.7 0.8 0.85 0.9 0.95 1.0];
    p.G2P  = linspace(0.5*p.P_amb, p.P_d*1.001, 500);
    p.G2P0 = p.P_amb + (p.P_d - p.P_amb)*fr;
    p.G2   = zeros(numel(fr), numel(p.G2P));
    for j = 1:numel(fr)
        [Pt, Gt, Gm] = build_hem(tab, p.G2P0(j), segs, p.ID1, K_ent, 0.5*p.P_amb, p.T_dl);
        p.G2(j,:) = interp1(Pt, Gt, p.G2P, 'linear', 0);  % 0 above the dewar pressure
        if j == numel(fr), p.in_Gmax = Gm; end
    end
    % Reverse path (tank liquid back to the dewar): 3/4 -> 1/2 -> valve -> 1/2 -> hose
    sgh  = (p.ID1/p.ID_hose)^2;
    if p.L3 > 0, Kc2 = 0.5*(1 - sg3); else, Kc2 = 0; end
    rsegs = struct('L',  {p.L3,     p.L2,     p.L1,                p.L_hose}, ...
                   'D',  {p.ID3,    p.ID2,    p.ID1,               p.ID_hose}, ...
                   'eps',{p.eps_ln, p.eps_ln, p.eps_ln,            p.eps_hose}, ...
                   'Kloc',{0,       0,        0,                   0}, ...
                   'Kin',{0,        Kc2,      p.K_valve + p.K_fit, (1-sgh)^2/sgh^2});
    if p.L_dw > 0
        K_dwd = 0;  if ~isnan(p.Cv_dewar), K_dwd = 894*(p.ID_dw/0.0254)^4/p.Cv_dewar^2; end
        rsegs(end+1) = struct('L',p.L_dw,'D',p.ID_dw,'eps',1.5e-6,'Kloc',0,'Kin',0.5*(1-sdw) + K_dwd);
    end
    rsegs = rsegs([rsegs.L] > 0);
    frr   = [0.5 0.7 0.85 0.95 1.0 1.1 1.25 1.5 2.0 3.0];
    p.GrP  = linspace(0.5*p.P_amb, (p.P_amb + (p.P_d - p.P_amb)*3)*1.001, 600);
    p.GrP0 = p.P_amb + (p.P_d - p.P_amb)*frr;
    p.Gr   = zeros(numel(frr), numel(p.GrP));
    for j = 1:numel(frr)
        [Pt, Gt] = build_hem(tab, p.GrP0(j), rsegs, p.ID1, 0.5, 0.5*p.P_amb);
        p.Gr(j,:) = interp1(Pt, Gt, p.GrP, 'linear', 0);
    end
    % Vapor vent: K_tot = f*L/D + 0.5 (entrance) + vent valve. f and gam are evaluated live in vent_flow,
    % so the Fanno table is built over a grid of gam and K_tot.
    if isnan(p.L_vent), p.L_vent = max(p.L_t - p.h_dip, 0.05); end
    if p.L_vent > 3
        warning('boiloff_v3:L_vent', ['p.L_vent = %.2f m. It is in metres: a length in inches needs *in2m ' ...
            '(a 0.152 in bore this long chokes the vent and stalls the fill).'], p.L_vent);
    end
    % vent valve: Cv on the tube velocity, plus reducer losses into / out of its bore if given
    K_vv = 894*(p.ID_vent/0.0254)^4/p.Cv_vent_valve^2;
    if ~isnan(p.ID_vent_valve)
        sv = min((p.ID_vent/p.ID_vent_valve)^2, 1);  K_vv = K_vv + (1-sv)^2 + 0.5*(1-sv);
    end
    p.K_vent_valve = K_vv;
    p.K_vent_minor = 0.5 + K_vv;
    p.f_vent0 = f_darcy(1e5, 1.5e-6/p.ID_vent);              % starting guess for the live friction factor
    p.vt_r   = linspace(0, 1, 401).';
    p.vt_gam = (1.25:0.025:1.75).';
    p.vt_lK  = log(logspace(log10(0.3), log10(2000), 60)).';
    PHI = zeros(numel(p.vt_r), numel(p.vt_gam), numel(p.vt_lK));
    for a = 1:numel(p.vt_gam)
        for b = 1:numel(p.vt_lK)
            [rr, ph] = build_vent(p.vt_gam(a), exp(p.vt_lK(b)));
            q = interp1(rr, ph, p.vt_r, 'linear');
            q(p.vt_r < rr(1)) = ph(1);  q(p.vt_r > rr(end)) = 0;   % choked below, no flow at r = 1
            PHI(:,a,b) = q;
        end
    end
    p.vt_F = griddedInterpolant({p.vt_r, p.vt_gam, p.vt_lK}, PHI, 'linear', 'nearest');
    % Two-phase vent: saturated liquid at P_tank through the vent tube to P_amb (HEM)
    if true
        vsegs = struct('L',p.L_vent, 'D',p.ID_vent, 'eps',1.5e-6, 'Kloc',0, 'Kin',0);
        P0s = linspace(1.03*p.P_amb, 8e5, 18);  Gl = zeros(size(P0s));
        for j = 1:numel(P0s)
            [Pt, Gt] = build_hem(tab, P0s(j), vsegs, p.ID_vent, 0.5 + p.K_vent_valve, 0.5*p.P_amb);
            Gl(j) = interp1(Pt, Gt, p.P_amb);
        end
        p.vl_P = [p.P_amb, P0s];  p.vl_G = [0, Gl];
    end
end

function f = f_darcy(Re, rel)
% Churchill friction factor: one smooth formula over laminar, transition and turbulent.
    Re = max(Re, 1);
    A = (2.457*log(1./((7./Re).^0.9 + 0.27*rel))).^16;
    B = (37530./Re).^16;
    f = 8*((8./Re).^12 + 1./(A + B).^1.5).^(1/12);
end

function [Ptab, Gtab, Gmax] = build_hem(tab, P0, segs, Dt, K_ent, Pmin, T_liq)
% Mass flux G (referenced to diameter Dt) of saturated liquid at P0 flashing through a series
% of segments, versus the outlet pressure (homogeneous equilibrium, adiabatic).
% For each trial G the pressure is marched along each segment (dz/dP from momentum + friction);
% stagnation enthalpy is conserved. Flow chokes where 1 + G^2 dv/dP -> 0.
% At each junction: Bernoulli for the area change plus the loss Kin on the new segment's velocity.
    nG = 80; nP = 400;
    T0 = Tsat_of_P(P0, tab);
    if nargin > 6 && ~isnan(T_liq), T0 = min(T_liq, T0); end   % subcooled source liquid
    r0 = sat_props(T0, tab);  h0 = r0.h_l;  vl0 = 1/r0.rho_l;
    Gs = logspace(log10(2), log10(4e4), nG);
    Pout = nan(1,nG);
    for i = 1:nG
        G  = Gs(i);
        G1 = G*(Dt/segs(1).D)^2;
        P  = P0 - (1 + K_ent)*G1^2*vl0/2;
        ok = (P > 1.05*Pmin);  vprev = vl0;  Gprev = G1;
        k = 1;
        while ok && k <= numel(segs)
            sg = segs(k);  Gk = G*(Dt/sg.D)^2;  D = sg.D;  L = sg.L;
            if k > 1
                P = P - ((1 + sg.Kin)*Gk^2 - Gprev^2)*vprev/2;
                if P <= 1.05*Pmin, ok = false; break, end
            end
            Pg = linspace(P, Pmin, nP).';
            rr = tab.F(Tsat_of_P(Pg, tab));
            rl = rr(:,2); rg = rr(:,3); hl = rr(:,6); hg = rr(:,7); mul = rr(:,13); mug = rr(:,14);
            vl = 1./rl; vg = 1./rg; vlg = vg - vl; hlg = hg - hl;
            a = Gk^2*vlg.^2/2;  b = Gk^2*vl.*vlg + hlg;  c = Gk^2*vl.^2/2 + hl - h0;
            x = (-b + sqrt(max(b.^2 - 4*a.*c, 0)))./(2*a);  x = min(max(x,0),1);
            v = vl + x.*vlg;  mum = 1./(x./mug + (1-x)./mul);  Re = Gk*D./mum;
            f = f_darcy(Re, sg.eps/D) + sg.Kloc*D/L;
            dvdP = gradient(v, Pg);
            term = 1 + Gk^2*dvdP;                             % -> 0 at choking
            integrand = term./(f*Gk^2.*v/(2*D));              % dz / d(-P)
            dz = 0.5*(integrand(2:end) + integrand(1:end-1)).*(Pg(1:end-1) - Pg(2:end));
            kc = find(term <= 0, 1);  if isempty(kc), kc = nP; end
            z  = [0; cumsum(dz)];  z = z(1:kc);
            if z(end) < L, ok = false; break, end
            [zu, iu] = unique(z);
            P = interp1(zu, Pg(iu), L);
            vprev = interp1(flipud(Pg), flipud(v), P);  Gprev = Gk;
            if P <= Pmin, ok = false; end
            k = k + 1;
        end
        if ok, Pout(i) = P; end
    end
    feas = ~isnan(Pout);
    if ~any(feas), error('boiloff_v3:flow','No feasible flow found in build_hem (P0 = %.0f Pa).', P0); end
    Gf = Gs(feas);  Po = Pout(feas);  Gmax = Gf(end);
    Ptab = [0.5*Pmin, fliplr(Po), P0];  Gtab = [Gmax, fliplr(Gf), 0];
    [Ptab, o] = sort(Ptab);  Gtab = Gtab(o);
    [Ptab, iu] = unique(Ptab);  Gtab = Gtab(iu);
end

function [r, phi] = build_vent(gam, Ktot)
% Fanno-flow vent table: Phi = G/sqrt(gam*P0*rho0) vs r = P_back/P0, friction K_tot.
    Ff  = @(M) (1-M.^2)./(gam*M.^2) + (gam+1)/(2*gam)*log((gam+1)*M.^2./(2+(gam-1)*M.^2));
    Mg  = linspace(1e-3, 0.99999, 4000);  Fg = Ff(Mg);
    PPs = @(M) (1./M).*sqrt((gam+1)./(2+(gam-1)*M.^2));
    M1c = interp1(fliplr(Fg), fliplr(Mg), Ktot);
    M1  = linspace(1e-3, M1c, 3000);
    F2  = Ff(M1) - Ktot;
    M2  = ones(size(M1));  ok = F2 > 1e-12;
    M2(ok) = interp1(fliplr(Fg), fliplr(Mg), F2(ok));
    P1P0 = (1+(gam-1)/2*M1.^2).^(-gam/(gam-1));
    r   = P1P0.*PPs(M2)./PPs(M1);
    phi = M1.*(1+(gam-1)/2*M1.^2).^(-(gam+1)/(2*(gam-1)));
    [r, o] = sort(r);  phi = phi(o);
    [r, iu] = unique(r);  phi = phi(iu);
end

function K = build_kappa(raw, fluid)
% Isentropic exponent of the gas on the gas table grid (T, P); saturated vapor at or below Tsat(P).
    Tg = raw.T;  Ps = raw.D(:,1);
    K = zeros(numel(raw.gasT), numel(raw.gasP));
    for j = 1:numel(raw.gasP)
        P = raw.gasP(j);
        if P < Ps(end), Ts = interp1(log(Ps), Tg, log(P), 'pchip'); else, Ts = Tg(end); end
        for i = 1:numel(raw.gasT)
            if raw.gasT(i) <= Ts + 0.5
                K(i,j) = cprops('ISENTROPIC_EXPANSION_COEFFICIENT','T',min(Ts,Tg(end)),'Q',1,fluid);
            else
                K(i,j) = cprops('ISENTROPIC_EXPANSION_COEFFICIENT','T',raw.gasT(i),'P',P,fluid);
            end
        end
    end
end

function M = mat_db(name)
% Wall material data. k and cp are NIST cryogenic-material curve fits:
%   log10(y) = sum a_n*(log10 T)^n,  n = 0..8   (trc.nist.gov/cryogenics/materials)
% k fits: 1-300 K; cp fits: 4-300 K. Both are held at their 300 K values above 300 K.
% rho and eps are typical handbook values (eps is an estimate; finish and frost change it).
    switch upper(strrep(strrep(name,' ',''),'-',''))
        case {'304SS','304','SS304','304STAINLESS'}
            M.name = '304 SS';  M.rho = 7900;  M.eps = 0.3;
            M.k  = [-1.4087 1.3982 0.2543 -0.6260 0.2334 0.4256 -0.4658 0.1650 -0.0199];
            M.cp = [22.0061 -127.5528 303.647 -381.0098 274.0328 -112.9212 24.7593 -2.239153 0];
        case {'6061AL','6061','AL6061','6061T6','6061T6AL'}
            M.name = '6061-T6 Al';  M.rho = 2700;  M.eps = 0.1;
            M.k  = [0.07918 1.0957 -0.07277 0.08084 0.02803 -0.09464 0.04179 -0.00571 0];
            M.cp = [46.6467 -314.292 866.662 -1298.3 1162.27 -637.795 210.351 -38.3094 2.96344];
        otherwise
            error('boiloff_v3:material','Unknown material ''%s''. Use ''304SS'' or ''6061Al''.', name);
    end
end

function k = k_mat(T, M)
% Wall thermal conductivity [W/m/K] from the material's NIST fit.
    k = 10.^polyval(fliplr(M.k), log10(min(max(T,1),300)));
end

function c = cp_mat(T, M)
% Wall specific heat [J/kg/K] from the material's NIST fit.
    c = 10.^polyval(fliplr(M.cp), log10(min(max(T,4),300)));
end

function C = cryo_db(name)
% Cryogen data: CoolProp fluid name, labels, property-cache file and the table temperature grids
% (saturation table from just above the triple point to a few K below the critical point).
    switch upper(name)
        case {'LOX','O2','OXYGEN'}
            C.name = 'LOX';  C.label = 'O2';  C.coolprop = 'Oxygen';  C.cache = 'lox_tables_v3.mat';
            C.Tsat_grid = 60:0.5:150;      % triple 54.4 K, critical 154.6 K
            C.Tvap_grid = 65:2.5:340;      % >= N2 triple point (pressurant table shares it)
        case {'LN2','N2','NITROGEN'}
            C.name = 'LN2';  C.label = 'N2';  C.coolprop = 'Nitrogen';  C.cache = 'ln2_tables_v3.mat';
            C.Tsat_grid = 64:0.5:122;      % triple 63.2 K, critical 126.2 K
            C.Tvap_grid = 64:2.5:339;
        otherwise
            error('boiloff_v3:cryogen','Unknown cryogen ''%s''. Use ''LOX'' or ''LN2''.', name);
    end
end

function T = Tsat_of_P(P, tab)
    P = min(max(P, tab.P(1)), tab.P(end));
    T = tab.TofLogP(log(P));
end

function s = sat_props(T, tab)
% Saturated cryogen properties at temperature T from the cached table.
    T = min(max(T, tab.T(1)), tab.T(end));
    r = tab.F(T);
    s.T=T; s.Psat=r(1); s.rho_l=r(2); s.rho_g=r(3); s.u_l=r(4); s.u_g=r(5);
    s.h_l=r(6); s.h_g=r(7); s.cp_l=r(8); s.cp_g=r(9); s.cv_g=r(10);
    s.k_l=r(11); s.k_g=r(12); s.mu_l=r(13); s.mu_g=r(14); s.sigma=r(15);
    s.beta_l=r(16); s.beta_g=r(17); s.kap=r(18);
end

function [k, mu, cp, rho, beta] = gas_film(Tf, P, p)
% Superheated cryogen gas properties at film temperature Tf and pressure P (table lookup).
    Tf = min(max(Tf, p.tab.gasT(1)), p.tab.gasT(end));
    P  = min(max(P,  p.tab.gasP(1)), p.tab.gasP(end));
    v  = p.tab.gas([Tf, log(P)]);  v = v(:);
    k = v(1); mu = v(2); cp = v(3); rho = v(4); beta = v(5);
end

function raw = build_tables(fluid, cryo, pgas)
    tb = tic;
% Saturated cryogen, superheated cryogen gas and air property tables from CoolProp (called once).
% Temperature ranges come from cryo_db (they must stay above the triple point).
    Tg = cryo.Tsat_grid(:);  n = numel(Tg);
    Dm = zeros(n,18);
    for i = 1:n
        T = Tg(i);
        Dm(i,:) = [ ...
            cprops('P','T',T,'Q',0,fluid), cprops('D','T',T,'Q',0,fluid), cprops('D','T',T,'Q',1,fluid), ...
            cprops('U','T',T,'Q',0,fluid), cprops('U','T',T,'Q',1,fluid), ...
            cprops('H','T',T,'Q',0,fluid), cprops('H','T',T,'Q',1,fluid), ...
            cprops('C','T',T,'Q',0,fluid), cprops('C','T',T,'Q',1,fluid), cprops('O','T',T,'Q',1,fluid), ...
            cprops('L','T',T,'Q',0,fluid), cprops('L','T',T,'Q',1,fluid), ...
            cprops('V','T',T,'Q',0,fluid), cprops('V','T',T,'Q',1,fluid), ...
            cprops('I','T',T,'Q',0,fluid), ...
            cprops('ISOBARIC_EXPANSION_COEFFICIENT','T',T,'Q',0,fluid), ...
            cprops('ISOBARIC_EXPANSION_COEFFICIENT','T',T,'Q',1,fluid), ...
            cprops('ISOTHERMAL_COMPRESSIBILITY','T',T,'Q',0,fluid) ];
    end
    raw.T = Tg;  raw.D = Dm;
    fprintf('  saturated %s table done (%s)\n', fluid, fmt_time(toc(tb)));

    % Superheated gas on (T, log P): [k mu cp rho beta]. Below saturation, saturated vapor.
    Tq = (60:5:330)';  Pq = exp(linspace(log(2e4), log(6e6), 40))';
    TsofP = @(P) interp1(log(Dm(:,1)), Tg, log(min(P, Dm(end,1))), 'pchip');
    props = {'L','V','C','D','ISOBARIC_EXPANSION_COEFFICIENT'};
    Gq = zeros(numel(Tq), numel(Pq), 5);
    for j = 1:numel(Pq)
        if mod(j, 10) == 0, fprintf('  %s gas table %d / %d (%s)\n', fluid, j, numel(Pq), fmt_time(toc(tb))); end
        if Pq(j) < Dm(end,1), Ts = TsofP(Pq(j)); else, Ts = Tg(end); end
        for i = 1:numel(Tq)
            for k = 1:5
                if Tq(i) <= Ts + 0.5
                    Gq(i,j,k) = cprops(props{k},'T',Ts,'Q',1,fluid);
                else
                    Gq(i,j,k) = cprops(props{k},'T',Tq(i),'P',Pq(j),fluid);
                end
            end
        end
    end
    raw.gasT = Tq;  raw.gasP = Pq;  raw.gas = Gq;
    raw.gasKap = build_kappa(raw, fluid);
    fprintf('  %s gas tables done (%s)\n', fluid, fmt_time(toc(tb)));

    Ta = (100:5:400)';  Da = zeros(numel(Ta),4);         % air at 1 atm: [k mu cp rho]
    for i = 1:numel(Ta)
        Da(i,:) = [ cprops('L','T',Ta(i),'P',101325,'Air'), cprops('V','T',Ta(i),'P',101325,'Air'), ...
                    cprops('C','T',Ta(i),'P',101325,'Air'), cprops('D','T',Ta(i),'P',101325,'Air') ];
    end
    raw.airT = Ta;  raw.airD = Da;

    % Internal energy and pressure on (T, rho) for the ullage state (two-phase inside the dome)
    Tv = cryo.Tvap_grid(:);  Rv = logspace(-2, log10(80), 80);
    Uv = zeros(numel(Tv), numel(Rv));  Pv = Uv;
    for i = 1:numel(Tv)
        for j = 1:numel(Rv)
            Uv(i,j) = cprops('U','T',Tv(i),'D',Rv(j),fluid);
            Pv(i,j) = cprops('P','T',Tv(i),'D',Rv(j),fluid);
        end
    end
    raw.vapT = Tv;  raw.vapR = Rv;  raw.vapU = Uv;  raw.vapP = Pv;
    fprintf('  %s ullage table done (%s)\n', fluid, fmt_time(toc(tb)));

    % ---- Pressurant: (T, rho) -> u, P on the same T grid (ullage mixture, bottles, fuel ullage) ----
    if strcmpi(pgas, 'Helium'), rmax = 150; else, rmax = 900; end
    Rp = logspace(-2, log10(rmax), 90);  Up = nan(numel(Tv), numel(Rp));  Pp = Up;
    for i = 1:numel(Tv)
        for j = 1:numel(Rp)
            try
                Up(i,j) = cprops('U','T',Tv(i),'D',Rp(j),pgas);  Pp(i,j) = cprops('P','T',Tv(i),'D',Rp(j),pgas);
            catch
            end
        end
    end
    Up = fillmissing(Up, 'nearest', 2);  Pp = fillmissing(Pp, 'nearest', 2);
    raw.pvR = Rp;  raw.pvU = Up;  raw.pvP = Pp;
    fprintf('  %s ullage table done (%s)\n', pgas, fmt_time(toc(tb)));
    % saturated pressurant liquid [Psat rho_l h_l] (dissolved-gas density and condensation enthalpy)
    if strcmpi(pgas, 'Helium')
        raw.psT = [2.2; 5];  raw.psD = [1 125 0; 2 125 0];
    else
        Tt = cprops1('Ttriple', pgas);  Tc = cprops1('Tcrit', pgas);
        Ts = (ceil(Tt)+1:0.5:floor(Tc)-4)';  Sd = zeros(numel(Ts),3);
        for i = 1:numel(Ts)
            Sd(i,:) = [cprops('P','T',Ts(i),'Q',0,pgas), cprops('D','T',Ts(i),'Q',0,pgas), cprops('H','T',Ts(i),'Q',0,pgas)];
        end
        raw.psT = Ts;  raw.psD = Sd;
    end
    % pressurant gas on (T, log P): [k mu cp rho beta h]; below saturation, saturated vapor
    Tq = (70:5:360)';  Pq = exp(linspace(log(2e4), log(2.5e7), 45))';
    props = {'L','V','C','D','ISOBARIC_EXPANSION_COEFFICIENT','H'};
    Gq = zeros(numel(Tq), numel(Pq), 6);
    for j = 1:numel(Pq)
        if mod(j, 10) == 0, fprintf('  %s gas table %d / %d (%s)\n', pgas, j, numel(Pq), fmt_time(toc(tb))); end
        Ts = -Inf;
        if ~strcmpi(pgas,'Helium') && Pq(j) < raw.psD(end,1)
            Ts = interp1(log(raw.psD(:,1)), raw.psT, log(Pq(j)), 'pchip');
        end
        for i = 1:numel(Tq)
            for k = 1:6
                if Tq(i) <= Ts + 0.5
                    Gq(i,j,k) = cprops(props{k},'T',Ts,'Q',1,pgas);
                else
                    Gq(i,j,k) = cprops(props{k},'T',Tq(i),'P',Pq(j),pgas);
                end
            end
        end
    end
    raw.pgT = Tq;  raw.pgP = Pq;  raw.pgD = Gq;
    fprintf('  %s gas table done (%s)\n', pgas, fmt_time(toc(tb)));
end

function [tab, air] = make_interpolants(raw)
% Build the interpolants once (much faster than interp1 'pchip', which refits every call).
    Dm = raw.D;
    tab.T = raw.T;  tab.P = Dm(:,1);  tab.rho_g = Dm(:,3);
    tab.rho_l_asc = flipud(Dm(:,2));
    tab.F        = griddedInterpolant(raw.T, Dm, 'pchip');                 % all 18 columns
    tab.Fres     = griddedInterpolant(raw.T, Dm(:,[2 3 4 5]), 'pchip');    % [rho_l rho_g u_l u_g]
    tab.TofLogP  = griddedInterpolant(log(Dm(:,1)), raw.T, 'pchip');
    tab.TofRhoG  = griddedInterpolant(Dm(:,3), raw.T, 'pchip');
    tab.TofRhoL  = griddedInterpolant(flipud(Dm(:,2)), flipud(raw.T), 'pchip');
    tab.ug       = griddedInterpolant(raw.T, Dm(:,5), 'pchip');
    tab.rhog     = griddedInterpolant(raw.T, Dm(:,3), 'pchip');
    tab.gas      = griddedInterpolant({raw.gasT, log(raw.gasP)}, raw.gas, 'linear', 'nearest');
    tab.gasT     = raw.gasT([1 end]);  tab.gasP = raw.gasP([1 end]);
    tab.gasMuKap = griddedInterpolant({raw.gasT, log(raw.gasP)}, cat(3, raw.gas(:,:,2), raw.gasKap), 'linear', 'nearest');
    tab.h_l_col  = Dm(:,6);                                              % for h_l -> T
    tab.vapT = raw.vapT;  tab.vlogR = log(raw.vapR(:));  tab.vapU = raw.vapU;  tab.vapP = raw.vapP;
    tab.plogR = log(raw.pvR(:));  tab.pvU = raw.pvU;  tab.pvP = raw.pvP;
    tab.psT = raw.psT([1 end]);   tab.psF = griddedInterpolant(raw.psT, raw.psD, 'pchip');
    tab.pgas = griddedInterpolant({raw.pgT, log(raw.pgP)}, raw.pgD, 'linear', 'nearest');
    tab.pgT = raw.pgT([1 end]);  tab.pgTv = raw.pgT;  tab.pgP = raw.pgP([1 end]);
    air.T = raw.airT;
    air.F = griddedInterpolant(raw.airT, raw.airD, 'pchip');
end

function v = cprops(out, in1, v1, in2, v2, fluid)
% Thin wrapper around CoolProp PropsSI through MATLAB's Python bridge.
    v = double(py.CoolProp.CoolProp.PropsSI(out, in1, v1, in2, v2, fluid));
end

function v = cprops1(out, fluid)
% CoolProp fluid constant (Ttriple, Tcrit, Pcrit, ...).
    v = double(py.CoolProp.CoolProp.PropsSI(out, fluid));
end
