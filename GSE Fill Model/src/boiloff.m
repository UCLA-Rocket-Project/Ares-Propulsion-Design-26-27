% Austin's Math
% Known Limitations: 
% - The tank is a thin vertical cylinder with end caps neglected
%   (The Mylar and insulation is also neglected).
% - The vent is modeled as adiabatic, line friction is not modeled and
%   pressure losses are lumped into one coefficient to match HF fill data.
% - Dewar does not have a CV in globe valve but model does not allow
%   backflow when tank pressure exceeds Dewar pressure (happens when lox
%   flashes at beginning of fill).
% - Venting Lox once plume is seen is not modeled, assumed that fill is
%   closed immediately and no liquid gets vented. 
% - Temperature of ullage and bulk liquid is assumed to be the same -> no
%   heat transfer between the two.  
% - Tank modeled with two axial nodes, wet (liquid) and dry (ullage). 
% - Dewar fill level and droop is not considered (should affect fill rate).
% - Aggressive boiling means we see plume before the liquid level actually 
%   reaches the vent. 
% - Condensation and ice forming on the outside of tanks and lines not
%   considered, affects natural convection from atmosphere to tank. 
% - Boiloff and temp. change during tank pressurization not considered. 
%
% The script outputs masses that are higher than we saw at hotfire, likely
% because of these limitations. However, it is fine for proof of concept,
% to show that additional vent cycles have diminishing returns. After 10
% it doesn't make much difference. Also proves that hitting the 25.6kg
% calculated at 90K is impossible. 

%% 1. Unit Conversion ----------------------------------------------------------
p = struct();           % start p empty every run (a script shares the base workspace)
psi2Pa = 6894.757;      % Pa per psi
in2m   = 0.0254;        % m per inch
ft2m   = 0.3048;        % m per foot

% ---- Tank (vertical cylinder, end caps neglected) ----
p.D        = 7.981*in2m;     % tank inner diameter                                 [m]
p.V_tank   = 1938.1*in2m^3;  % total internal tank volume                          [m^3]
p.t_w      = 0.322*in2m;     % tank wall thickness                                 [m]
p.V_dip    = 1370.4*in2m^3;  % liquid volume when the level reaches the dip-tube [m^3]
p.k_w      = 12;             % 304 SS conductivity, for conduction between the wet and dry wall [W/m/K]

% ---- Vent (dip tube) ----
p.vtube_OD   = 0.25*in2m;    % vent / dip tube outer diameter (-4)                   [m]
p.vtube_wall = 0.049*in2m;   % vent tube wall                                        [m]
p.Cd_vent    = 0.66;         % vent discharge coefficient (adjusted to fit HF data)

% ---- Fill path: dewar -> hose -> tube 1 -> ball valve -> tube 2 -> tube 3 (3/4 in) -> tank bottom ----
% I don't remember the actual length of these tubes so they are guesses. 
p.L_hose   = 5*ft2m;         % insulated cryogenic flex hose (-8 AN), treated as adiabatic [m]
p.ID_hose  = 0.5*in2m;       % Assumed hose ID;                                        [m]
p.eps_hose = 1.5e-6;         % Assumed hose roughness                                  [m]
p.OD1 = 0.5*in2m;  p.t1 = 0.065*in2m;  p.L1 = 1.5*ft2m;   % uninsulated 1/2 in tube, dewar side of the valve
p.OD2 = 0.5*in2m;  p.t2 = 0.065*in2m;  p.L2 = 1.0*ft2m;   % uninsulated 1/2 in tube, valve -> 3/4 in tube
p.OD3 = 0.75*in2m; p.t3 = 0.065*in2m;  p.L3 = 1.5*ft2m;   % uninsulated 3/4 in tube -> tank bottom
p.eps_ln     = 1.5e-6;       % Assumed tube roughness                                  [m]
p.K_ent      = 0.5;          % entrance loss from the dewar
p.Cv_valve   = 9.7;          % fill ball valve Cv, fully open
p.valve_open = 0.9;          % fraction open (0-1)
p.valve_curve = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1; ...      % opening fraction
                 0 .02 .05 .09 .15 .23 .34 .48 .65 .83 1]; % Cv fraction: typical full-port ball valve 
                                                            % (couldn't find official from valworx)                                                   
p.Cv_dewar   = 1.1;          % dewar liquid globe valve from HF, fully open. RegO lists
                             % Cv 1.10 for its 3/8 in valves, 0.70 for 1/4 in. or use NaN = not modeled.
% Vacuum jacketed (assumed), treated as adiabatic. 0 = not modeled.
p.L_dw  = 0;               % Assumed diptube length     [m]
p.ID_dw = 0.305*in2m;  % Assumed dewar dip tube size 3/8 in OD -> 0.305 ID       [m]
p.K_fit      = 0;            % extra lumped loss (AN fittings, tee, bends), on the 1/2 in tube velocity

% ---- Operating conditions ----
p.P_d_psig = 110;            % dewar pressure, held constant for the whole fill      [psig]
p.T_amb    = 286;            % ambient temperature                                   [K]
p.P_amb    = 92389.75;       % ambient pressure                                      [Pa]
p.P_0      = NaN;            % initial tank pressure; NaN = P_amb                    [Pa]
p.T_wall0  = NaN;            % initial tank wall temperature; NaN = T_amb            [K]
p.T_ln1_0  = NaN;            % initial temperature of tube 1 (dewar side of the valve).
                             %   NaN = Tsat(P_dewar): it sits full of dewar liquid for a few
                             %   minutes before the fill valve is opened.               [K]

% ---- Procedure ----
p.n_fill   = 10;              % number of fill phases (A), with a vent-only phase (B) between them.
                             % Assumes no liquid ever gets vented (not true)
                    
p.t_vent   = 120;            % duration of each vent-only phase B
p.tA_max   = 3600;           % cap on one fill phase                          [s]

p.make_plots = true;         
p.rel_tol    = 1e-3;         % ODE relative tolerance


% ---- Heat transfer constants ----
p.q_ext_mult = 1.0;          % multiplier on outside heat leak (raise for frost)
% Boiling, wall -> liquid:  q = h*dT,  dT = T_wall - Tsat
%   h = h_film above the Leidenfrost superheat dT_leid (vapor film on the wall), h_nb below it.
%   dT_leid rises with pressure: dT_leid = dT_leid0 + dT_leid_slope*(P[bar] - 3.5).
%   Values fitted using HF2 and HF4 data
p.h_film        = 180;      % film boiling coefficient                              [W/m^2/K]
p.h_nb          = 1300;     % effective nucleate / transition boiling coefficient    [W/m^2/K]
p.dT_leid0      = 89;       % Leidenfrost superheat at 3.5 bar abs                  [K]
p.dT_leid_slope = 18.5;     % rise of the Leidenfrost superheat with pressure        [K/bar]
p.eps       = 0.3;           % emissivity of the outside of 304 SS (Condensation on tanks
                                % changes this significantly. Not modeled).
p.rho_w     = 7900;          % 304 SS density                                        [kg/m^3]
p.gam_v     = 1.49;          % O2 vapor gamma for the vent orifice 

p.fluid = 'Oxygen';

%% 2. DERIVED CONSTANTS AND PROPERTY TABLES --------------------------------
here = fileparts(mfilename('fullpath'));
p.g     = 9.80665;                       % gravity [m/s^2]
p.sigSB = 5.670374e-8;                   % Stefan-Boltzmann [W/m^2/K^4]

p.A_c   = pi*p.D^2/4;                    % tank cross-section
p.L_t   = p.V_tank/p.A_c;                % equivalent cylinder length
p.h_dip = p.V_dip/p.A_c;                 % dip-tube inlet height
p.H_wet = min(p.h_dip, p.L_t);           % wet wall node: tank bottom up to the dip-tube inlet
p.H_dry = max(p.L_t - p.H_wet, 1e-3);    % dry wall node: dip-tube inlet to the top
p.A_wet = pi*p.D*p.H_wet;   p.m_wet = p.rho_w*p.t_w*p.A_wet;   % inner area [m^2], mass [kg]
p.A_dry = pi*p.D*p.H_dry;   p.m_dry = p.rho_w*p.t_w*p.A_dry;
p.G_ax  = p.k_w*pi*p.D*p.t_w/(0.5*p.L_t);    % wet <-> dry conductance [W/K]

p.ID_vent = p.vtube_OD - 2*p.vtube_wall; 
p.A_v     = pi*p.ID_vent^2/4;
for s = '123'
    p.(['ID' s])   = p.(['OD' s]) - 2*p.(['t' s]);
    p.(['m_ln' s]) = p.rho_w*pi/4*(p.(['OD' s])^2 - p.(['ID' s])^2)*p.(['L' s]);
end
p.A_in   = pi*p.ID1^2/4;                 % inflow G is referenced to the 1/2 in tube ID
Cv_eff   = p.Cv_valve*interp1(p.valve_curve(1,:), p.valve_curve(2,:), p.valve_open);
p.K_valve = 894*(p.ID1/in2m)^4/Cv_eff^2; % K = 894 d^4/Cv^2, tuned to match HF fill data

if isnan(p.P_0), p.P_0 = p.P_amb; end
if isnan(p.T_wall0), p.T_wall0 = p.T_amb; end
p.P_d = p.P_d_psig*psi2Pa + p.P_amb;     % dewar absolute pressure [Pa]

% Property tables (cached next to this file, rebuilt if stale)
TAB_VERSION = 4;                          % 4: saturated table trimmed to the 9 columns used
cacheFile = fullfile(here, 'lox_tables_simple.mat');
tabRaw = [];
if exist(cacheFile, 'file')
    S = load(cacheFile);
    if isfield(S,'tabRaw') && isfield(S.tabRaw,'version') && S.tabRaw.version == TAB_VERSION ...
            && strcmp(S.tabRaw.fluid, p.fluid)
        tabRaw = S.tabRaw;
    else
        fprintf('Property cache is out of date; rebuilding.\n');
    end
end
if isempty(tabRaw)
    fprintf('Building property tables with CoolProp (one time, ~1 min)...\n');
    tabRaw = build_tables(p.fluid);
    tabRaw.version = TAB_VERSION;  tabRaw.fluid = p.fluid;
    save(cacheFile, 'tabRaw');
end
[p.tab, p.air] = make_interpolants(tabRaw);

if isnan(p.T_ln1_0), p.T_ln1_0 = Tsat_of_P(p.P_d, p.tab); end
p.Tsat_d = Tsat_of_P(p.P_d, p.tab);
p.sd     = sat_props(p.Tsat_d, p.tab);

% Flow tables
p = setup_flows(p);
fprintf('Fill valve: Cv %.1f at %.0f%% open -> Cv %.2f, K = %.2f\n', p.Cv_valve, 100*p.valve_open, Cv_eff, p.K_valve);
fprintf('Inflow: choked G = %.0f kg/m2/s at full dewar pressure (max %.3f kg/s)\n', p.in_Gmax, p.in_Gmax*p.A_in);
r_crit = (2/(p.gam_v+1))^(p.gam_v/(p.gam_v-1));
fprintf('Vent: tube ID %.3f in, Cd %.2f; choked above %.2f bar abs\n', p.ID_vent/in2m, p.Cd_vent, p.P_amb/r_crit/1e5);
fprintf('Wall: wet node %.1f kg (%.2f m), dry node %.1f kg (%.2f m)\n', p.m_wet, p.H_wet, p.m_dry, p.H_dry);

%% 3. INITIAL CONDITIONS ---------------------------------------------------
rho0 = cprops('D','T',p.T_amb,'P',p.P_0,p.fluid);
M0   = rho0*p.V_tank;                                       % initial O2 mass (gas)
U0   = M0*cprops('U','T',p.T_amb,'P',p.P_0,p.fluid);
% x = [M; U; T_ln1; T_ln2; T_ln3; m_vent; m_in; T_wet; T_dry]
x0 = [M0; U0; p.T_ln1_0; p.T_amb; p.T_amb; 0; 0; p.T_wall0; p.T_wall0];

%% 4. PROCEDURE ------------------------------------------------------------
absTol = [1e-4; 100; 1e-2; 1e-2; 1e-2; 1e-4; 1e-4; 1e-2; 1e-2];
opts   = odeset('RelTol',p.rel_tol,'AbsTol',absTol,'MaxStep',20);
optsA  = odeset(opts,'Events',@(t,x) ev_fill_full(t,x,p));

t = 0;  x = x0;
T_all = t;  X_all = x.';  ph_all = 1;          % 1 = A (fill), 2 = B (vent)
plog = cell(0,8);  N_fill_done = 0;
tv = p.t_vent(:);
for k = 1:p.n_fill
    % ---- Phase A: fill + vent until the liquid reaches the dip-tube inlet ----
    [tA,xA,~,~,ieA] = ode15s(@(tt,xx) rhs(tt,xx,p,true), [t t+p.tA_max], x, optsA);
    [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tA,xA,1);
    t = tA(end);  x = xA(end,:).';  N_fill_done = k;
    plog(end+1,:) = phase_line(sprintf('fill %d', k), t, x, p); 
    if isempty(ieA)
        warning('boiloff:fill', ['Fill %d did not reach the dip tube within tA_max (%.0f s). ' ...
            'Check P_dewar vs tank pressure, vent size and V_dip.'], k, p.tA_max);
        break
    end
    if k == p.n_fill, break, end                % last fill ends the run
    % ---- Phase B: fill closed, vent open for t_vent ----
    [tB,xB] = ode15s(@(tt,xx) rhs(tt,xx,p,false), [t t+tv(min(k,numel(tv)))], x, opts);
    [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tB,xB,2);
    t = tB(end);  x = xB(end,:).';
    plog(end+1,:) = phase_line(sprintf('vent %d', k), t, x, p); 
end

%% 5. POST-PROCESSING ------------------------------------------------------
n = numel(T_all);
D = struct('T',zeros(n,1),'P',zeros(n,1),'Vl',zeros(n,1), ...
           'ml',zeros(n,1),'mdot_in',zeros(n,1),'mdot_vent',zeros(n,1), ...
           'Qwl',zeros(n,1),'Qwg',zeros(n,1));
for i = 1:n
    [~,d] = rhs(T_all(i), X_all(i,:).', p, ph_all(i)==1);
    D.T(i)=d.T;  D.P(i)=d.P;  D.Vl(i)=d.Vl;  D.ml(i)=d.ml;
    D.mdot_in(i)=d.mdot_in;  D.mdot_vent(i)=d.mdot_vent;
    D.Qwl(i)=d.Qwl;  D.Qwg(i)=d.Qwg;
end
stEnd = fluid_state(X_all(end,1), X_all(end,2), p);
Tw_end = X_all(end,8:9);                   % [wet dry]

fprintf('\n%-11s %8s %10s %9s %9s %9s %9s %9s\n','phase end','t [s]','Tsat [K]','P [psig]','M [kg]','liq [kg]','in [kg]','vent [kg]');
for i = 1:size(plog,1)
    fprintf('%-11s %8.0f %10.1f %9.1f %9.2f %9.2f %9.2f %9.2f\n', plog{i,:});
end
fprintf('\n================ RESULT ================\n');
fprintf('Fill phases                    : %d\n', N_fill_done);
fprintf('Vent-open time (fill / vent)   : %.0f s / %.0f s\n', sum(diff(T_all).*(ph_all(2:end)==1)), ...
    sum(diff(T_all).*(ph_all(2:end)>1)));
Tsat_end = Tsat_of_P(stEnd.P, p.tab);
fprintf('Final tank pressure            : %.1f psig  (%.2f bar abs)\n', (stEnd.P - p.P_amb)/psi2Pa, stEnd.P/1e5);
fprintf('Saturation temp at that P      : %.2f K   <- Tsat(P_tank)\n', Tsat_end);
if stEnd.two_phase
    fprintf('Final liquid temperature       : %.2f K   (equilibrium model: liquid = vapor = Tsat)\n', stEnd.T);
else
    fprintf('Final gas temperature          : %.2f K   (no liquid in the tank)\n', stEnd.T);
end
fprintf('Final liquid density           : %.1f kg/m^3\n', stEnd.rho_l);
fprintf('O2 in tank when vent closes    : %.2f kg  (liquid %.2f kg)\n', X_all(end,1), stEnd.ml);
fprintf('O2 loaded from dewar           : %.2f kg\n', X_all(end,7));
fprintf('O2 vented (vapor only)         : %.2f kg\n', X_all(end,6));
fprintf('Dewar pressure (held)          : %.1f psig\n', p.P_d_psig);
fprintf('Peak tank pressure             : %.1f psig at t = %.0f s\n', (max(D.P) - p.P_amb)/psi2Pa, T_all(find(D.P == max(D.P),1)));
fprintf('Mass balance residual          : %.2e kg (should be ~0)\n', X_all(end,1) - M0 - X_all(end,7) + X_all(end,6));
fprintf('Final wall temps (wet / dry)   : %.1f K / %.1f K\n', Tw_end(1), Tw_end(2));
fprintf('========================================\n');

% ---- Plots ----
if p.make_plots
    tm = T_all/60;
    figure('Name','LOX chill-in','Position',[80 80 1300 780]);
    subplot(2,3,1);
    plot(tm, D.T, 'LineWidth',1.5); hold on;
    plot(tm, X_all(:,8:9), 'LineWidth',1.1);
    xlabel('time [min]'); ylabel('T [K]'); grid on;
    legend({'fluid (= T_{sat})','wall, wet (below dip tube)','wall, dry (above)'},'Location','best');
    subplot(2,3,2);
    plot(tm, (D.P - p.P_amb)/psi2Pa, 'LineWidth',1.5, 'DisplayName','tank'); hold on;
    yline(p.P_d_psig, '--k', 'DisplayName','dewar');
    xlabel('time [min]'); ylabel('P [psig]'); grid on; legend('Location','best');
    subplot(2,3,3);
    plot(tm, D.Vl/p.V_dip, 'LineWidth',1.5); hold on;
    yline(1,'--r','dip tube');
    xlabel('time [min]'); ylabel('V_{liquid} / V_{dip}'); grid on;
    subplot(2,3,4);
    plot(tm, X_all(:,1), tm, X_all(:,6), tm, X_all(:,7), 'LineWidth',1.5);
    xlabel('time [min]'); ylabel('mass [kg]'); grid on; legend('in tank','vented','loaded','Location','best');
    subplot(2,3,5);
    plot(tm, D.mdot_in, tm, D.mdot_vent, 'LineWidth',1.2);
    xlabel('time [min]'); ylabel('mass flow [kg/s]'); grid on; legend('in','vent (vapor)','Location','best');
    subplot(2,3,6);
    plot(tm, D.Qwl/1e3, tm, D.Qwg/1e3, 'LineWidth',1.2);
    xlabel('time [min]'); ylabel('wall heat [kW]'); grid on; legend('wall -> liquid','wall -> gas','Location','best');
end

%% ========================= FUNCTIONS ================================

function [T_all,X_all,ph_all] = append_hist(T_all,X_all,ph_all,tt,xx,code)
    T_all  = [T_all; tt(2:end)];  X_all = [X_all; xx(2:end,:)];
    ph_all = [ph_all; code*ones(numel(tt)-1,1)];
end

function c = phase_line(name, t, x, p)
    st = fluid_state(x(1), x(2), p);
    c = {name, t, Tsat_of_P(st.P, p.tab), (st.P - p.P_amb)/6894.757, x(1), st.ml, x(7), x(6)};
end

function [dx, d] = rhs(~, x, p, fill)
% RHS  Time derivative of the state vector and diagnostics.
%   x = [M; U; T_ln1; T_ln2; T_ln3; m_vent; m_in (loaded); T_wet; T_dry]
%   fill = true in phase A (fill valve open), false otherwise. The vent is always open.
    M = x(1); U = x(2); T1 = x(3); T2 = x(4); T3 = x(5);
    Twet = x(8);  Tdry = x(9);
    g = p.g;  sp = @(z) 0.5*(z + sqrt(z^2 + 1e-8));            % smooth max(0,z)

    % ---- Tank fluid state from (M, U) ----
    st = fluid_state(M, U, p);  T = st.T;  P = st.P;
    if st.two_phase, s = st.s; Tgas = T; else, s = []; Tgas = st.Tg; end

    % ---- Wall: wet node (partly covered by liquid) and dry node (vapor only) ----
    hl  = min(max(st.Vl/p.A_c, 0), p.L_t);               % liquid level
    fw  = min(hl/p.H_wet, 1);                              % covered fraction of the wet node
    Lg  = max(p.L_t - hl, 1e-2);                           % height of wall in the ullage
    A_wg = (1-fw)*p.A_wet;                                 % uncovered part of the wet node
    Tws  = (A_wg*Twet + p.A_dry*Tdry)/(A_wg + p.A_dry);    % mean dry-surface temperature
    [kg,mug,cpg,rhog,bg] = gas_film(0.5*(Tws + Tgas), P, p);   % ullage gas at film temperature
    Qwet_l = 0;
    if st.two_phase && fw > 0
        Qwet_l = fw*p.A_wet*boil_q(Twet - T, s, p);
    end
    Qwet_g = A_wg   *nat_conv(kg,mug,cpg,rhog,bg,Twet-Tgas,Lg,g)*(Twet - Tgas);
    Qdry_g = p.A_dry*nat_conv(kg,mug,cpg,rhog,bg,Tdry-Tgas,Lg,g)*(Tdry - Tgas);
    Qex_wet = p.A_wet*ext_heat(Twet, p, p.L_t, 'v');
    Qex_dry = p.A_dry*ext_heat(Tdry, p, p.L_t, 'v');
    Qax = p.G_ax*(Tdry - Twet);                            % conduction dry -> wet through the wall
    Qwl = Qwet_l;  Qwg = Qwet_g + Qdry_g;

    % ---- Fill line ----
    Qex_1 = pi*p.OD1*p.L1*ext_heat(T1, p, p.OD1, 'h');
    Qex_2 = pi*p.OD2*p.L2*ext_heat(T2, p, p.OD2, 'h');
    Qex_3 = pi*p.OD3*p.L3*ext_heat(T3, p, p.OD3, 'h');
    Pd = p.P_d;  sd = p.sd;                                % dewar held at P_d, saturated liquid
    % Tube 1 is always full of dewar liquid at Tsat(P_dewar).
    Q1 = pi*p.ID1*p.L1*boil_q(T1 - sd.T, sd, p);
    % Tubes 2 and 3 see the flashing stream (fill) or tank fluid (valve closed).
    Tloc = Tsat_of_P(P, p.tab);  sl = sat_props(Tloc, p.tab);
    Q2 = seg_heat(T2, p.ID2, p.L2, st.two_phase || fill, sl, Tloc, Tgas, P, p);
    Q3 = seg_heat(T3, p.ID3, p.L3, st.two_phase || fill, sl, Tloc, Tgas, P, p);

    % ---- Fill-line flow (dewar -> tank only; no flow when the tank is at or above dewar pressure) ----
    mdot_in = 0;  md0 = 0;
    if fill
        dP0 = 2000;                                           % linear below dP = 0.3 psi (flow ~ sqrt(dP) has infinite slope at 0)
        dP  = Pd - P;
        if dP > 0
            Pb  = min(max(Pd - max(dP,dP0), p.G_P(1)), p.G_P(end));
            md0 = p.A_in*interp1(p.G_P, p.G_G, Pb)*min(dP, dP0)/dP0;
        end
        % Line heat adds vapor and cuts the forward flow. For friction-dominated flow md^2*vbar
        % is fixed; heat spread along the line raises the mean specific volume by dx/2*v_lg:
        %   md^2 + bq*md - md0^2 = 0,  bq = Q*v_lg/(2*vbar*h_fg)
        Pm  = 0.5*(Pd + P);  sm = sat_props(Tsat_of_P(Pm, p.tab), p.tab);
        vl  = 1/sm.rho_l;  vlg = 1/sm.rho_g - vl;  hfg_m = sm.h_g - sm.h_l;
        x_ex = min(max((sd.h_l - sl.h_l)/(sl.h_g - sl.h_l), 0), 1);   % flash quality at the tank
        vbar = vl + 0.5*x_ex*vlg;
        bq  = (sp(Q1) + sp(Q2) + sp(Q3))*vlg/(2*vbar*hfg_m);
        mdot_in = 0.5*(-bq + sqrt(bq^2 + 4*md0^2));
        % tube-1 heat goes with the stream when it flows, otherwise back to the dewar
        Hdot_in = mdot_in*sd.h_l + Q1*mdot_in/(mdot_in + 0.005) + Q2 + Q3;
    else
        Hdot_in = Q2 + Q3;    % valve closed: tube 2/3 heat soaks into the tank fluid
    end

    % ---- Vent ----
    if st.two_phase
        hv = s.h_g;  rho0 = s.rho_g;
    else
        % Single-phase vapor: look up by (T, density), never (T, P), because a state on the
        % saturation line makes (T, P) ambiguous and CoolProp errors.
        Tv = st.Tg;  rho0 = vapor_side_rho(Tv, st.rho, p);
        try
            hv = cprops('H','T',Tv,'D',rho0,p.fluid);
        catch
            sv = sat_props(Tv, p.tab);  hv = sv.h_g;
        end
    end
    mvent = vent_flow(P, rho0, p);                        % vapor only: the level never reaches the vent

    % ---- Conservation ----
    dM = mdot_in - mvent;
    dU = Qwl + Qwg + Hdot_in - mvent*hv;

    % ---- Wall and line temperatures ----
    dTwet = (Qex_wet - Qwet_l - Qwet_g + Qax)/(p.m_wet*cw304(Twet));
    dTdry = (Qex_dry - Qdry_g - Qax)/(p.m_dry*cw304(Tdry));
    dT1 = (Qex_1 - Q1)/(p.m_ln1*cw304(T1));
    dT2 = (Qex_2 - Q2)/(p.m_ln2*cw304(T2));
    if p.L3 > 0, dT3 = (Qex_3 - Q3)/(p.m_ln3*cw304(T3)); else, dT3 = 0; end

    dx = [dM; dU; dT1; dT2; dT3; mvent; mdot_in; dTwet; dTdry];
    d = struct('T',T,'P',P,'Vl',st.Vl,'ml',st.ml,'mdot_in',mdot_in,'mdot_vent',mvent, ...
               'Qwl',Qwl,'Qwg',Qwg);
end

function Q = seg_heat(Tw, ID, L, wetted, sl, Tloc, Tgas, P, p)
% Heat from an uninsulated fill-line segment into the fluid inside it.
    if L <= 0, Q = 0; return; end
    A = pi*ID*L;
    if wetted
        Q = A*boil_q(Tw - Tloc, sl, p);
    else
        [k,mu,cp,rho,b] = gas_film(0.5*(Tw + Tgas), P, p);
        Q = A*nat_conv(k,mu,cp,rho,b,Tw-Tgas,ID,p.g)*(Tw - Tgas);
    end
end

function v = Vliq(x, p)
% Liquid volume in the tank.
    st = fluid_state(x(1), x(2), p);  v = st.Vl;
end

function [value, isterminal, direction] = ev_fill_full(~, x, p)
% Liquid reaches the dip-tube inlet: end of a fill phase.
    value = Vliq(x, p) - p.V_dip;  isterminal = 1;  direction = 1;
end

function st = fluid_state(M, U, p)
% FLUID_STATE  Tank state from total mass M and internal energy U.
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
            error('boiloff:state','Single-phase state lookup failed (rho = %.3g kg/m3, u = %.4g J/kg): %s', ...
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

function q = boil_q(dTe, s, p)
% Wall -> liquid heat flux [W/m^2] at wall superheat dTe = T_wall - T_sat:
%   q = h*dTe, h = h_film above the Leidenfrost superheat, h_nb below it
%   (switch smoothed over ~3 K so the ODE solver sees no jump).
    dTL = max(p.dT_leid0 + p.dT_leid_slope*(s.Psat/1e5 - 3.5), 10);
    h   = p.h_nb + (p.h_film - p.h_nb)/(1 + exp(-(dTe - dTL)/3));
    q   = h*dTe;
end

function q = ext_heat(Tw, p, L, kind)
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
    q = p.q_ext_mult*(h*dT + p.eps*p.sigSB*(p.T_amb^4 - Tw^4));
end

function mdot = vent_flow(P, rho0, p)
% Vapor vent as a compressible orifice (isentropic nozzle with discharge coefficient Cd):
%   mdot = Cd*A*sqrt(gam*P*rho0)*phi(r),  r = P_amb/P
%   choked  (r <= r_crit): phi = (2/(gam+1))^((gam+1)/(2(gam-1)))
%   subsonic           : phi = sqrt(2/(gam-1)*(r^(2/gam) - r^((gam+1)/gam)))
    dP  = P - p.P_amb;
    dPe = 0.5*(dP + sqrt(dP^2 + 5^2));                % smooth max(dP, 0)
    r   = p.P_amb/(p.P_amb + dPe);
    gm  = p.gam_v;  rc = (2/(gm+1))^(gm/(gm-1));
    if r <= rc
        phi = (2/(gm+1))^((gm+1)/(2*(gm-1)));
    else
        phi = sqrt(2/(gm-1)*(r^(2/gm) - r^((gm+1)/gm)));
    end
    mdot = p.Cd_vent*p.A_v*phi*sqrt(gm*P*rho0);
end

function p = setup_flows(p)
% Build the inflow table G(P_tank) at the fixed dewar pressure.
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
                  'Kin',{0,       K_hose,     Kc + p.K_fit, p.K_valve, (1-sg3)^2/sg3^2});   % Borda-Carnot expansion
    segs = segs([segs.L] > 0);
    % One table at the (fixed) dewar pressure: mass flux G vs tank pressure
    [Pt, Gt, p.in_Gmax] = build_hem(tab, p.P_d, segs, p.ID1, K_ent, 0.5*p.P_amb);
    p.G_P = linspace(0.5*p.P_amb, p.P_d, 500);
    p.G_G = interp1(Pt, Gt, p.G_P, 'linear', 0);          % 0 at the dewar pressure
end

function f = f_darcy(Re, rel)
% Churchill friction factor: one smooth formula over laminar, transition and turbulent.
    Re = max(Re, 1);
    A = (2.457*log(1./((7./Re).^0.9 + 0.27*rel))).^16;
    B = (37530./Re).^16;
    f = 8*((8./Re).^12 + 1./(A + B).^1.5).^(1/12);
end

function [Ptab, Gtab, Gmax] = build_hem(tab, P0, segs, Dt, K_ent, Pmin)
% Mass flux G (referenced to diameter Dt) of saturated liquid at P0 flashing through a series
% of segments, versus the outlet pressure (homogeneous equilibrium, adiabatic).
% For each trial G the pressure is marched along each segment (dz/dP from momentum + friction);
% stagnation enthalpy is conserved. Flow chokes where 1 + G^2 dv/dP -> 0.
% At each junction: Bernoulli for the area change plus the loss Kin on the new segment's velocity.
    nG = 80; nP = 400;
    T0 = Tsat_of_P(P0, tab);  r0 = sat_props(T0, tab);  h0 = r0.h_l;  vl0 = 1/r0.rho_l;
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
            rl = rr(:,2); rg = rr(:,3); hl = rr(:,6); hg = rr(:,7); mul = rr(:,8); mug = rr(:,9);
            vl = 1./rl; vg = 1./rg; vlg = vg - vl; hlg = hg - hl;
            a = Gk^2*vlg.^2/2;  b = Gk^2*vl.*vlg + hlg;  c = Gk^2*vl.^2/2 + hl - h0;
            x = (-b + sqrt(max(b.^2 - 4*a.*c, 0)))./(2*a);  x = min(max(x,0),1);
            v = vl + x.*vlg;  mum = 1./(x./mug + (1-x)./mul);  Re = Gk*D./mum;
            f = f_darcy(Re, sg.eps/D);
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
    if ~any(feas), error('boiloff:flow','No feasible flow found in build_hem (P0 = %.0f Pa).', P0); end
    Gf = Gs(feas);  Po = Pout(feas);  Gmax = Gf(end);
    Ptab = [0.5*Pmin, fliplr(Po), P0];  Gtab = [Gmax, fliplr(Gf), 0];
    [Ptab, o] = sort(Ptab);  Gtab = Gtab(o);
    [Ptab, iu] = unique(Ptab);  Gtab = Gtab(iu);
end

function c = cw304(T)
% 304 SS specific heat [J/kg/K], NIST cryogenic-materials fit (4-300 K, ~5%); held above 300 K.
    a = [22.0061 -127.5528 303.647 -381.0098 274.0328 -112.9212 24.7593 -2.239153];
    x = log10(min(max(T,4),300));
    c = 10.^polyval(fliplr(a), x);
end

function T = Tsat_of_P(P, tab)
    P = min(max(P, tab.P(1)), tab.P(end));
    T = tab.TofLogP(log(P));
end

function s = sat_props(T, tab)
% Saturated O2 properties at temperature T from the cached table.
    T = min(max(T, tab.T(1)), tab.T(end));
    r = tab.F(T);
    s.T=T; s.Psat=r(1); s.rho_l=r(2); s.rho_g=r(3); s.h_l=r(6); s.h_g=r(7);
end

function rho = vapor_side_rho(T, rho_in, p)
% Keep a single-phase vapor density just on the vapor side of saturation.
    rho = rho_in;
    if T <= p.tab.T(end)
        rho = min(rho_in, 0.999*p.tab.rhog(T));
    end
end

function [k, mu, cp, rho, beta] = gas_film(Tf, P, p)
% Superheated O2 gas properties at film temperature Tf and pressure P (table lookup).
    Tf = min(max(Tf, p.tab.gasT(1)), p.tab.gasT(end));
    P  = min(max(P,  p.tab.gasP(1)), p.tab.gasP(end));
    v  = p.tab.gas([Tf, log(P)]);  v = v(:);
    k = v(1); mu = v(2); cp = v(3); rho = v(4); beta = v(5);
end

function raw = build_tables(fluid)
% Saturated O2, superheated O2 gas and air property tables from CoolProp (called once).
    Tg = (60:0.5:150)';  n = numel(Tg);
    % Saturated columns: [P rho_l rho_g u_l u_g h_l h_g mu_l mu_g]
    Dm = zeros(n,9);
    for i = 1:n
        T = Tg(i);
        Dm(i,:) = [ ...
            cprops('P','T',T,'Q',0,fluid), cprops('D','T',T,'Q',0,fluid), cprops('D','T',T,'Q',1,fluid), ...
            cprops('U','T',T,'Q',0,fluid), cprops('U','T',T,'Q',1,fluid), ...
            cprops('H','T',T,'Q',0,fluid), cprops('H','T',T,'Q',1,fluid), ...
            cprops('V','T',T,'Q',0,fluid), cprops('V','T',T,'Q',1,fluid) ];
    end
    raw.T = Tg;  raw.D = Dm;

    % Superheated gas on (T, log P): [k mu cp rho beta]. Below saturation, saturated vapor.
    Tq = (60:5:330)';  Pq = exp(linspace(log(2e4), log(6e6), 40))';
    TsofP = @(P) interp1(log(Dm(:,1)), Tg, log(min(P, Dm(end,1))), 'pchip');
    props = {'L','V','C','D','ISOBARIC_EXPANSION_COEFFICIENT'};
    Gq = zeros(numel(Tq), numel(Pq), 5);
    for j = 1:numel(Pq)
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

    Ta = (100:5:400)';  Da = zeros(numel(Ta),4);         % air at 1 atm: [k mu cp rho]
    for i = 1:numel(Ta)
        Da(i,:) = [ cprops('L','T',Ta(i),'P',101325,'Air'), cprops('V','T',Ta(i),'P',101325,'Air'), ...
                    cprops('C','T',Ta(i),'P',101325,'Air'), cprops('D','T',Ta(i),'P',101325,'Air') ];
    end
    raw.airT = Ta;  raw.airD = Da;
end

function [tab, air] = make_interpolants(raw)
% Build the interpolants once (much faster than interp1 'pchip', which refits every call).
    Dm = raw.D;
    tab.T = raw.T;  tab.P = Dm(:,1);  tab.rho_g = Dm(:,3);
    tab.rho_l_asc = flipud(Dm(:,2));
    tab.F        = griddedInterpolant(raw.T, Dm, 'pchip');                 % all 9 columns
    tab.Fres     = griddedInterpolant(raw.T, Dm(:,[2 3 4 5]), 'pchip');    % [rho_l rho_g u_l u_g]
    tab.TofLogP  = griddedInterpolant(log(Dm(:,1)), raw.T, 'pchip');
    tab.TofRhoG  = griddedInterpolant(Dm(:,3), raw.T, 'pchip');
    tab.TofRhoL  = griddedInterpolant(flipud(Dm(:,2)), flipud(raw.T), 'pchip');
    tab.ug       = griddedInterpolant(raw.T, Dm(:,5), 'pchip');
    tab.rhog     = griddedInterpolant(raw.T, Dm(:,3), 'pchip');
    tab.gas      = griddedInterpolant({raw.gasT, log(raw.gasP)}, raw.gas, 'linear', 'nearest');
    tab.gasT     = raw.gasT([1 end]);  tab.gasP = raw.gasP([1 end]);
    air.T = raw.airT;
    air.F = griddedInterpolant(raw.airT, raw.airD, 'pchip');
end

function v = cprops(out, in1, v1, in2, v2, fluid)
% Thin wrapper around CoolProp PropsSI through MATLAB's Python bridge.
    v = double(py.CoolProp.CoolProp.PropsSI(out, in1, v1, in2, v2, fluid));
end
