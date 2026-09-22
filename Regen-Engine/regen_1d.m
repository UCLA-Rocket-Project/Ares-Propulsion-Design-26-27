% Test script for 1d regen engine

clc; clear; close all;
set(groot, 'defaultFigureColor', 'w', 'defaultAxesColor', 'w', ...
    'defaultAxesXColor', 'k', 'defaultAxesYColor', 'k', 'defaultTextColor', 'k');

%% Basic Parameters
F = 1465 * 4.44822; % Target thrust in N
Pamb = 12.3; % psia
Param.T_amb = 298; % K (from feed calc sheet)
mfrac_eth = 0.75; mfrac_h2o = 1 - mfrac_eth; % fuel ratio
mw_eth = 46.068; mw_h2o = 18.015; % g/mol
x_eth = (mfrac_eth/mw_eth) / (mfrac_eth/mw_eth + mfrac_h2o/mw_h2o); x_h2o = 1 - x_eth; % Mole fraction
P_crit = [6140000, 22064000]; % ethanol, water
Param.P_crit_mix = x_eth * P_crit(1) + x_h2o * P_crit(2); % fuel mixture
oxidizer = 'LOX';
% Assumed efficiencies
cstar_eff = 0.90;
cf_eff = 0.98;

% CEA parameters
o_f = 1.1;
Pc_us = 375; % psia, target
Param.Pc = convpres(Pc_us, 'psi', 'Pa');
card_str = sprintf(['fuel C2H5OH(L)   C 2 H 6 O 1\n', ...
    'h,cal=-66370.0      t(k)=298.00      wt%%=75.00\n', ...
    'fuel water H 2.0 O 1.0  wt%%=25.00\n', ...
    'h,cal=-68308.  t(k)=298.00 rho,g/cc = 0.9998']);

py.rocketcea.cea_obj.add_new_fuel('ETHANOL_WATER_75_25(L)', card_str);
fuel = 'ETHANOL_WATER_75_25(L)';
c = py.rocketcea.cea_obj.CEA_Obj(pyargs('oxName', oxidizer,'fuelName', fuel));
exp_ratio = c.get_eps_at_PcOvPe(pyargs('Pc', Pc_us, 'MR', o_f, 'PcOvPe',(Pc_us / Pamb)));

% Cantera and read yaml
ct = py.importlib.import_module('cantera');
Cantera = ct.Solution('nasa9_species.yaml');
Y_fuel_toal = 1/(1+ o_f);
Y_ox_total = o_f / (1 + o_f);
Y_eth = 0.75 * Y_fuel_toal;
Y_water = 0.25 * Y_fuel_toal;
Y_o2 = Y_ox_total;
Y_str = sprintf('C2H5OH: %.8f, H2O: %.8f, O2:%.8f', Y_eth, Y_water, Y_o2);

% transport properties [Cp, mu, k, Prandtl]
transport_chamber = c.get_Chamber_Transport(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
transport_throat = c.get_Throat_Transport(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
transport_exit = c.get_Exit_Transport(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
i_stag = c.get_Enthalpies(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
vel_sonic = c.get_SonicVelocities(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));

cp_g_ref = [double(transport_chamber{1}), double(transport_throat{1}), double(transport_exit{1})] .* 4184; % J/kg*K
mu_g_ref = [double(transport_chamber{2}), double(transport_throat{2}), double(transport_exit{2})] .* 0.0001; % Pa*s viscosity
k_g_ref = [double(transport_chamber{3}), double(transport_throat{3}), double(transport_exit{3})] .* 0.4184; % W/m*K thermal conductivity
prandtl_ref = [double(transport_chamber{4}), double(transport_throat{4}), double(transport_exit{4})];
Param.cp_g = cp_g_ref(1); Param.mu_g = mu_g_ref(1); Param.k_g = k_g_ref(1); Param.prandtl = prandtl_ref(1);
i_stag_ref = [double(i_stag{1}), double(i_stag{2}), double(i_stag{3})] .* 2326; % J/kg specific enthalpy
vel_sonic_ref = [double(vel_sonic{1}), double(vel_sonic{2}), double(vel_sonic{3})] .* 0.3048; % m/s local speed of sound

gamma_chamber = c.get_Chamber_MolWt_gamma(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
Gas.gamma = double(gamma_chamber{2}); % specific heat ratio, constant for isentropic relations

temps_cea = c.get_Temperatures(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio)); % Mat.r
Gas.T_stag = double(temps_cea{1}) * 5/9 * cstar_eff^2; % K, Chamber (adiabatic Wall -> constant stagnation temp)

cstar_theo = double(c.get_Cstar(pyargs('Pc', Pc_us, 'MR', o_f))) * 0.3048; % m/s
cf_cea = c.get_PambCf(pyargs('Pamb', Pamb, 'Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
cf_theo = double(cf_cea{1});
Param.MW = 23.446; % g/mol
Param.MW_coolant = x_eth*mw_eth + x_h2o*mw_h2o; % g/mol, coolant blend

exit_pres_ratio = c.get_PcOvPe(pyargs('Pc', Pc_us, 'MR', o_f, 'eps', exp_ratio));
Gas.P_exit = convpres(Pc_us/exit_pres_ratio, 'psi', 'Pa');

% Calculating target injector manifold pressure
fuel_stiffness = 0.20; % standard
dP_inj = fuel_stiffness * Param.Pc; % Pa
mdot_total = F / (cstar_theo * cstar_eff * cf_theo * cf_eff); % kg/s
Param.cstar_act = cstar_theo .* cstar_eff; % m/s
Param.mdot_f = mdot_total ./ (1 + o_f); % kg/s
rho_f_inj = 789 * mfrac_eth + 1000 * mfrac_h2o; % kg/m^3 (20C)
CdA_f_inj = Param.mdot_f / sqrt(2 * rho_f_inj * dP_inj); % m^2 (Heritage)
P_target = Param.Pc + dP_inj; % Pa manifold pressure (target value)

%% Throat Geometry (Rao)
Geo.At = (Param.cstar_act * mdot_total)/Param.Pc; % m^2
Geo.Rt = sqrt(Geo.At/pi); % m
Geo.R_curve = 1.5 * Geo.Rt;

%% Diverging Geometry
Geo.Ae = Geo.At * exp_ratio; % m^2
Geo.R_exit = sqrt(Geo.Ae/pi); % m
Geo.percent_len = 0.8; % input (0.8 is optimal fractional length for most cases)
Geo.len_div = Geo.percent_len * ((Geo.R_exit - Geo.Rt)/tan(deg2rad(15))); % m diverging length from radii
Geo.theta_e = deg2rad(13); Geo.theta_n = deg2rad(23); % deg, from HH fig 4.16

%% Chamber and Converging Geometry
Geo.conv_angle = deg2rad(45); % deg, standard
Geo.L_star = 35; % in, optimal for ethanol/lox
Geo.CR = 5; % A_c / A_t
Geo.V_total = convlength(Geo.L_star, 'in', 'm') * Geo.At; % m^3
Geo.Rc = Geo.Rt * sqrt(Geo.CR);
Geo.x_conv_tangent = -(Geo.R_curve) * sin(Geo.conv_angle);
Geo.y_conv_tangent = Geo.Rt + (Geo.R_curve) * (1 - cos(Geo.conv_angle));
% 1.5*Geo.Rt converging fillet
Geo.R_fillet = 1.5*Geo.Rt;
Geo.y_conv_fillet = Geo.Rc - Geo.R_fillet; % center height of fillet circle
Geo.y_f_tangent = Geo.y_conv_fillet + Geo.R_fillet*cos(Geo.conv_angle);
Geo.x_f_tangent = Geo.x_conv_tangent - (Geo.y_f_tangent - Geo.y_conv_tangent) / tan(Geo.conv_angle);
Geo.x_conv_fillet = Geo.x_f_tangent - Geo.R_fillet * sin(Geo.conv_angle); % x center of fillet circle
Geo.x_chamber_end = Geo.x_conv_fillet;

Geo.len_conv = -Geo.x_chamber_end; % m
Geo.V_conv = (1/3 * pi * Geo.len_conv) * (Geo.Rc^2 + Geo.Rt^2 + Geo.Rc*Geo.Rt); % m^3
Geo.V_chamber = Geo.V_total - Geo.V_conv;
Geo.len_chamber = Geo.V_chamber/(pi*Geo.Rc^2);
Geo.len_total = Geo.len_chamber + Geo.len_conv + Geo.len_div; % m

%% Encoding Geometry
Geo.dx = 0.001; % 1 mm step size
Geo.x_t = 0; % Recentering coordinates to around throat
Geo.x_chamber_start = -(Geo.len_chamber + Geo.len_conv);
Geo.x_exit = Geo.len_div;
% Position arrrays
Geo.pos_i = Geo.x_chamber_start:Geo.dx:Geo.x_exit; % axial position
Geo.pos_j = zeros(size(Geo.pos_i)); % radii
Geo.pos_throat = floor((Geo.len_chamber + Geo.len_conv)/Geo.dx);
Geo.pos_conv = floor((Geo.len_chamber)/Geo.dx);
Geo.pos_taper = 0; % Use if tapering geometry towards the injector
% Diverging arc
Geo.xn = (0.382 * Geo.Rt) * sin(Geo.theta_n); % diverging tangent point where parabola starts (Geo.theta_n)
Geo.Rn = Geo.Rt + (0.382 * Geo.Rt) * (1 - cos(Geo.theta_n)); % y coordinate of Geo.xn
% Rao parabola (y = ax^2 + bx + c)
% Need to derive parabola that hits (Geo.xn, Geo.Rn), (xe, Geo.Re) w/ starting slope tan(Geo.theta_n)
Geo.matrix_A = [Geo.xn^2, Geo.xn, 1; Geo.x_exit^2, Geo.x_exit, 1; 2*Geo.xn, 1, 0];
Geo.matrix_B = [Geo.Rn; Geo.R_exit; tan(Geo.theta_n)];
Geo.coeffs = Geo.matrix_A \ Geo.matrix_B;
Geo.a_Rao = Geo.coeffs(1); b_Rao = Geo.coeffs(2); c_Rao = Geo.coeffs(3);

% Encoder loop
x = Geo.pos_i;
m1 = x < Geo.x_chamber_end;
m2 = x >= Geo.x_chamber_end & x < Geo.x_f_tangent;
m3 = x >= Geo.x_f_tangent & x < Geo.x_conv_tangent;
m4 = x >= Geo.x_conv_tangent & x < 0;
m5 = x >= 0 & x < Geo.xn;
m6 = x >= Geo.xn;
Geo.pos_j(m1) = Geo.Rc;
Geo.pos_j(m2) = Geo.y_conv_fillet + sqrt(Geo.R_fillet^2 - (x(m2) - Geo.x_conv_fillet).^2);
Geo.pos_j(m3) = Geo.y_f_tangent - tan(Geo.conv_angle) .* (x(m3) - Geo.x_f_tangent);
Geo.pos_j(m4) = (Geo.Rt + Geo.R_curve) - sqrt((Geo.R_curve)^2 - x(m4).^2);
Geo.pos_j(m5) = (Geo.Rt + 0.382*Geo.Rt) - sqrt((0.382*Geo.Rt)^2 - x(m5).^2);
Geo.pos_j(m6) = Geo.a_Rao .* x(m6).^2 + b_Rao .* x(m6) + c_Rao;

% Slant Length
Geo.dx_slope = gradient(Geo.pos_j, Geo.dx);
Geo.dl = Geo.dx .* sqrt(1 + Geo.dx_slope .^2);

%% Cooling Channel and Outer Jacket Geometry
Geo.min_tol = 0.001; % m, 3d printer minimum feature
Geo.D_gas = 2 * Geo.pos_j; % Array of gas-side diameter Geo.At every node
Geo.D_t = 2 * Geo.Rt;
Geo.coat_thickness = 0.0002; % m, Thermal Coating Thickness
Geo.out_wall_thickness = 0.001; % m, Outer Jacket
Geo.w_rib = Geo.min_tol; % fixed, rib width
Geo.min_cw = Geo.min_tol;

% Variable Channel Height
Geo.h_channel_chamber_forward = Geo.min_tol*10;
Geo.h_channel_chamber_aft = Geo.min_tol*7;
Geo.h_channel_throat = Geo.min_tol*1.5;
Geo.h_channel_nozzle = Geo.min_tol*3;
% Variable Inner Wall Thickness
Geo.wall_thickness_forward = Geo.min_tol*1.25;
Geo.wall_thickness_aft = Geo.min_tol*1;
Geo.wall_thickness_throat = Geo.min_tol*1;
Geo.wall_thickness_nozzle = Geo.min_tol*2;

m_conv = 1:Geo.pos_conv;
m_thr = (Geo.pos_conv + 1):Geo.pos_throat;
m_div = (Geo.pos_throat + 1):length(Geo.pos_i);
% Interpolate Wall Thickness
Geo.wall_thickness = zeros(size(Geo.pos_i));
Geo.wall_thickness(m_conv) = -((Geo.wall_thickness_forward - Geo.wall_thickness_aft)/Geo.pos_conv).*(m_conv-1) + Geo.wall_thickness_forward;
Geo.wall_thickness(m_thr)  = ((Geo.wall_thickness_aft - Geo.wall_thickness_throat)/(Geo.D_gas(Geo.pos_conv) - Geo.D_gas(Geo.pos_throat))) .* ...
    (Geo.D_gas(m_thr) - Geo.D_gas(Geo.pos_conv)) + Geo.wall_thickness_aft;
Geo.wall_thickness(m_div)  = ((Geo.wall_thickness_nozzle - Geo.wall_thickness_throat)/(Geo.D_gas(end) - Geo.D_gas(Geo.pos_throat))) .* ...
    (Geo.D_gas(m_div) - Geo.D_gas(Geo.pos_throat)) + Geo.wall_thickness_throat;

Geo.D_channel_base = Geo.D_gas + 2 * Geo.wall_thickness; % Engine diameters with added wall thickness
Geo.D_channel_base_throat = Geo.D_channel_base(Geo.pos_throat);
Geo.D_channel_base_chamber = Geo.D_channel_base(1);
Geo.D_channel_base_nozzle = Geo.D_channel_base(length(Geo.pos_i));
Geo.D_t_base = Geo.D_t + 2 * Geo.wall_thickness_throat; % Throat diameter with added wall thickness
% Number of channels determined Geo.At throat
Geo.circ_t_base = pi * (Geo.D_t_base); % Gas side diameter + wall thickness
Geo.num_channel = floor(Geo.circ_t_base / (Geo.w_rib + Geo.min_cw));
Geo.circ_local_base = pi * Geo.D_channel_base; % Local circumferences across engine
Geo.w_channel = (Geo.circ_local_base - (Geo.num_channel * Geo.w_rib)) ./ Geo.num_channel; % variable, channel widths

% Interpolate channel height
Geo.h_channel = zeros(size(Geo.pos_i));
Geo.h_channel(m_conv) = -((Geo.h_channel_chamber_forward - Geo.h_channel_chamber_aft)/Geo.pos_conv).*(m_conv-1) + Geo.h_channel_chamber_forward;
Geo.h_channel(m_thr)  = ((Geo.h_channel_chamber_aft - Geo.h_channel_throat)/(Geo.D_channel_base_chamber - Geo.D_channel_base_throat)) .* ...
    (Geo.D_channel_base(m_thr) - Geo.D_channel_base_chamber) + Geo.h_channel_chamber_aft;
Geo.h_channel(m_div)  = ((Geo.h_channel_throat - Geo.h_channel_nozzle)/(Geo.D_channel_base_throat - Geo.D_channel_base_nozzle)) .* ...
    (Geo.D_channel_base(m_div) - Geo.D_channel_base_throat) + Geo.h_channel_throat;

Geo.D_h = (4 .* Geo.w_channel .* Geo.h_channel) ./ (2 * Geo.w_channel + 2 * Geo.h_channel); % hydraulic diameter of rectangular channels
Geo.Per_heated = Geo.w_channel + 2 * Geo.h_channel; % heated perimeter
% HT Areas
Geo.A_gas = pi .* Geo.D_gas .* Geo.dl; % Gas-wall convection SA
Geo.A_w = pi .* ((Geo.D_gas + Geo.D_channel_base)./2) .* Geo.dl; % wall-wall conduction SA (average diameter)
Geo.A_co = Geo.Per_heated .* Geo.num_channel .* Geo.dl; % coolant side surface area (heated)
Geo.A_wc = Geo.w_channel .* Geo.dl; % cool wall area per increment

%% Visualization Plot % CAD Geometry Export
r_channel_base = Geo.D_channel_base ./ 2;
r_outer_jacket = r_channel_base + Geo.h_channel;
r_outer_wall = r_channel_base + Geo.h_channel + Geo.out_wall_thickness;

figure('Name', '1D Engine Geometry', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Geo.pos_j, 'k', 'LineWidth', 2, 'DisplayName', 'Hot Wall');
plot(Geo.pos_i, -Geo.pos_j, 'k', 'LineWidth', 2, 'HandleVisibility','off');
plot(Geo.pos_i, r_channel_base, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Cold Wall');
plot(Geo.pos_i, -r_channel_base, 'b--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(Geo.pos_i, r_outer_jacket, 'b', 'LineWidth', 2, 'DisplayName', 'Outer Jacket');
plot(Geo.pos_i, -r_outer_jacket, 'b', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(Geo.pos_i, r_outer_wall, 'k', 'LineWidth', 2, 'DisplayName', 'Outer Wall');
plot(Geo.pos_i, -r_outer_wall, 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
title('Regen 1D Profile')
xlabel('Axial Position x (m)');
ylabel('Radial Position y (m)');
axis equal;
xline(0, 'r--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'geometry.pdf', 'ContentType','vector');
hold off;

zero_array = zeros(length(Geo.pos_i), 1);
writematrix([zero_array, r_channel_base', Geo.pos_i'], 'channel_base.txt');
writematrix([zero_array, r_outer_jacket', Geo.pos_i'], 'outer_jacket.txt');
writematrix([zero_array, r_outer_wall', Geo.pos_i'], 'outer_wall.txt');
writematrix([zero_array, Geo.pos_j', Geo.pos_i'], 'hot_wall.txt');

%% Gas and Coolant Properties
% 1D interpolation from 3 CEA points for Cp, Gas.gamma, k, mu for Bartz
Gas.M_local = zeros(size(Geo.pos_i));
Gas.A_local = Geo.D_gas.^2 .* (pi/4);
Gas.AR_local = Gas.A_local ./ Geo.At;

for k = 1:length(Geo.pos_i)
    if Geo.pos_i(k) < 0 % Chamber and Converging
        Gas.M_local(k) = flowisentropic(Gas.gamma, Gas.AR_local(k), 'sub');
    elseif Geo.pos_i(k) == 0
        Gas.M_local(k) = 1.0;
    else
        Gas.M_local(k) = flowisentropic(Gas.gamma, Gas.AR_local(k), 'sup');
    end
end
Gas.M_chamber =  flowisentropic(Gas.gamma, Geo.Rc^2/Geo.Rt^2, 'sub');
Gas.M_ref = [Gas.M_chamber, 1.0, Gas.M_local(end)];

Gas.pressure = Param.Pc * (1 + (Gas.gamma - 1)/2 .* Gas.M_local .^2) .^ (-Gas.gamma / (Gas.gamma - 1)); % From local Mach relations
Gas.vel_sonic_local = interp1(Gas.M_ref, vel_sonic_ref, Gas.M_local, 'linear', 'extrap');
Gas.velocity_local = Gas.vel_sonic_local .* Gas.M_local;
Gas.cp_g_local = interp1(Gas.M_ref, cp_g_ref, Gas.M_local, 'linear', 'extrap');
Gas.mu_g_local = interp1(Gas.M_ref, mu_g_ref, Gas.M_local, 'linear', 'extrap');
Gas.k_g_local = interp1(Gas.M_ref, k_g_ref, Gas.M_local, 'linear', 'extrap');
Gas.prandtl_g_local = interp1(Gas.M_ref, prandtl_ref, Gas.M_local, 'linear', 'extrap');
Gas.i_stag_local = interp1(Gas.M_ref, i_stag_ref, Gas.M_local, 'linear', 'extrap');

Gas.Taw = Gas.T_stag * ((1 + Gas.prandtl_g_local.^(1/3).*((Gas.gamma - 1) / 2) .* Gas.M_local.^2) ...
    ./ (1 + ((Gas.gamma - 1) / 2) .* Gas.M_local.^2));
Gas.i_aw = Gas.i_stag_local .* (1 + (Gas.prandtl_g_local.^(1/3) - 1) .* (Gas.velocity_local.^2 ./ 2) ./ Gas.i_stag_local);

Cool = load_coolprop('coolprop_tables.mat', 3, x_eth, x_h2o, mfrac_eth, mfrac_h2o);

%% Material Properties
Mat.r = 1e-4; % surface roughness

Mat.k_w_ref_temps = [22 233 448 657 866 1079 1289 1500] + 273.15;
Mat.k_w_ref = [11.9 13.7 16.9 21.7 25.6 22.9 19.1 17.7];
Mat.ref_temps = [21 537 815 982 1093] + 273.15;
Mat.E_ref = [16.5 15.2 11 5.5 3.4] .* 10^10;
Mat.nu_ref = [0.3 0.28 0.323 0.368 0.4];
Mat.alpha_ref_temps = [427 538 649 871 1093 1827] + 273.15;
Mat.alpha_ref = [1.436 1.49 1.543 1.745 1.834 1.834] .* 10*-5;

% Thermal Coating: Type YSZ Zirconium Oxide
Mat.k_Al2O3_ref_temps = [300 400 500 600 700 800 900 1000 1100 1200 1300 ...
    1400 1500 1600];
Mat.k_Al2O3_ref = [36.9601 27.209 20.93 16.3045 12.558 10.4650 8.9999 7.9534 7.1792 ...
    6.8023 6.2467 5.9651 5.8604 5.8604];
Mat.k_ZrO2_ref_temps = [100 1400] + 273.15;
Mat.k_ZrO2_ref = [2.93 2.512];

%% Output Arrays
arr_names = {'T_fin', 'T_ow', 'T_hw', 'T_tc', 'T_cw', 'T_bulk', 'P', 'P_loss', ...
             'q_flux', 'h_i', 'h_c', 'h_c_f', 'h_tp', 'h_nb', 'CHF', 'CHF_safe', ...
             'stress', 'q_eq', 'f', 'Nu', 'T_sat', 'Re', 'vel_c', 'fin_eff', ...
             'q_error', 'iter_T', 'sigma', 'h_c_FEA', 'h_g', 'gas_flux', 'k_coating'};
for idx = 1:length(arr_names)
    Arrays.([arr_names{idx}, '_array']) = zeros(size(Geo.pos_i));
end

%% Main Loop
Loop.P_guess = convpres(500, 'psi', 'Pa');
Loop.P_prev_guess = convpres(450, 'psi', 'Pa');
Loop.P_error = realmax;
Loop.P_prev_error = 0;
Loop.tol_P = 100; % Pa
Loop.iter_P = 0;

while abs(Loop.P_error) > Loop.tol_P % Pressure guess loop
    Loop.iter_P = Loop.iter_P + 1;
    Loop.P_loc = Loop.P_guess;
    Loop.P_reduced = Loop.P_loc / Param.P_crit_mix;
    Loop.T_bulk = Param.T_amb;

    for d = length(Geo.pos_i):-1:1 % Axial marching loop
        % Local Geometry, add channel height array earlier
        disp(d);
        Loop.cw = Geo.w_channel(d);
        Loop.ch = Geo.h_channel(d);
        Loop.A_conv = Geo.num_channel * (Loop.cw + 2*Loop.ch) * Geo.dl(d);
        Loop.A_g_loc = Geo.A_gas(d);
        Loop.A_w_loc = Geo.A_w(d);
        Loop.A_base_loc = pi * Geo.D_channel_base(d) * Geo.dl(d);
        Loop.D_g_loc = Geo.D_gas(d);
        Loop.D_h_loc = (2*Loop.cw*Loop.ch)/(Loop.cw+Loop.ch); % hydraulic diameter
        Loop.A_c_cs = Loop.cw*Loop.ch; % Cross sectional area of channel
        Loop.A_c_cs_tot = Loop.A_c_cs * Geo.num_channel;
        if (d ~= 1) % If not Geo.At final chamber station
            Loop.A_c_cs_next = Geo.w_channel(d-1)*Loop.ch;
            Loop.D_h_loc_next = (2*Geo.w_channel(d-1)*Loop.ch/(Geo.w_channel(d-1)+Loop.ch));
        else % Geo.At final station
            Loop.A_c_cs_next = Loop.A_c_cs;
            Loop.D_h_loc_next = Loop.D_h_loc;
        end

        % Coolant properties
        Loop.rho_c = Cool.get_rho(Loop.P_loc, Loop.T_bulk);
        Loop.cp_c = Cool.get_cp(Loop.P_loc, Loop.T_bulk);
        Loop.mu_c = Cool.get_mu(Loop.P_loc, Loop.T_bulk);
        Loop.k_c = Cool.get_k(Loop.P_loc, Loop.T_bulk);
        Loop.T_sat = Cool.get_T_sat(Loop.P_loc);
        Loop.T_sat_200 = Cool.get_T_sat(convpres(200, 'psi', 'Pa'));
        Loop.rho_c_l = Cool.get_rho_l(Loop.P_loc);
        Loop.rho_c_v = Cool.get_rho_v(Loop.P_loc);
        Loop.surften = Cool.get_surften(Loop.P_loc);
        Loop.h_fg = Cool.get_h_fg(Loop.P_loc);
        Loop.quality = max(0,(Loop.cp_c * (Loop.T_bulk - Loop.T_sat)) / Loop.h_fg); % capped at 0 for Elfaham
        Loop.prandtl_c = Loop.mu_c * Loop.cp_c / Loop.k_c;
        Loop.vel_c = Param.mdot_f / (Loop.A_c_cs_tot * Loop.rho_c);
        Loop.Re = (Loop.rho_c * Loop.D_h_loc * Loop.vel_c) / Loop.mu_c;
        Arrays.Re_array(d) = Loop.Re;
        Arrays.vel_c_array(d) = Loop.vel_c;

        if (Loop.Re>4000) % Turbulent
            f_guess = 64/Loop.Re; % Initial guess
            f_error = realmax;
            tol_f = 0.0001;
            iter_F = 0;
            while abs(f_error) > tol_f
                iter_F = iter_F + 1;
                Loop.f_calc = 1/(-2*log10(2.5226/(Loop.Re*sqrt(f_guess))+Mat.r/(Loop.D_h_loc*3.7065)))^2;
                f_error = Loop.f_calc - f_guess;
                f_guess = 0.5*f_guess + 0.5*Loop.f_calc; % avoid overshooting
                if iter_F > 100
                    break;
                end
            end
        elseif (Loop.Re<2300) % Laminar
            Loop.f_calc = 64/Loop.Re;
        else % Transition
            Loop.f_calc = 0.02783 + (7.17*10^-6)*(Loop.Re - 2300);
        end
        Arrays.f_array(d) = Loop.f_calc;
        
        Loop.Pr = (Loop.cp_c*Loop.mu_c)/Loop.k_c;
        Loop.Nu = ((Loop.f_calc/8)*(Loop.Re-1000)*Loop.Pr)/...
            (1+12.7*(Loop.f_calc/8)^0.5*(Loop.Pr^(2/3)-1)); % Gnielinski 
        Loop.h_c = Loop.Nu * Loop.k_c / Loop.D_h_loc; % Coolant convection htc hydraulic diameter
        Arrays.h_c_array(d) = Loop.h_c;
        Arrays.Nu_array(d) = Loop.Nu;
        
        % Gas properties
        Loop.Taw_loc = Gas.Taw(d); % Local adiabatic wall temp
        
        %Temp Loops
        Temp = temp_iteration(Param, Cantera, Y_str, Geo, Gas, Cool, Mat, Loop, d);
        Temp = outer_wall_temp(Geo, Loop, Temp, Param, d);
        
        Arrays.T_fin_array(d) = Temp.T_fin_tip;
        Arrays.T_ow_array(d) = Temp.T_ow;
        Arrays.T_cw_array(d) = Temp.T_cw;
        Arrays.T_tc_array(d) = Temp.T_tc_guess;
        Arrays.T_hw_array(d) = Temp.T_hw;
        Arrays.fin_eff_array(d) = Temp.fin_eff;
        Arrays.h_i_array(d) = Temp.h_i;
        Arrays.q_error_array(d) = Temp.q_error;
        Arrays.iter_T_array(d) = Temp.iter_T;
        Arrays.q_eq_array(d) = Temp.q_eq;
        Arrays.h_c_f_array(d) = Temp.h_c_f;
        Arrays.h_nb_array(d) = Temp.h_nb;
        Arrays.h_tp_array(d) = Temp.h_tp;
        Arrays.sigma_array(d) = Temp.sigma;
        % FEA Arrays
        Arrays.h_c_FEA_array(d) = Temp.h_c_FEA;
        Arrays.h_g_array(d) = Temp.h_g;
        Arrays.gas_flux_array(d) = Temp.q_flux_gas;
        Arrays.k_coating_array(d) = Temp.k_coating;
        
        % Check CHF
        q_flux = Temp.q_eq/Loop.A_conv;
        Arrays.q_flux_array(d) = q_flux;
        dT_sub = Loop.T_sat - Loop.T_bulk; % K, bulk subcooling
        F_p = 1.17-8.56*(10^(-4))*convpres(Loop.P_loc, 'Pa', 'psi');
        CHF_base = 0.1003+0.05264*sqrt(convvel(Loop.vel_c, 'm/s', 'ft/s')*(dT_sub * 9/5)); 
        Arrays.CHF_array(d) = CHF_base*F_p*1635000;
        Arrays.CHF_array_safe(d) = Arrays.CHF_array(d)*0.9;
        Arrays.T_sat_array(d) = Loop.T_sat;

        % Prepare for next station
        Arrays.T_bulk_array(d) = Loop.T_bulk; % Store the updated bulk temperature
        Loop.T_bulk = Loop.T_bulk + Temp.q_eq/(Param.mdot_f*Loop.cp_c); % K

        % Calculate pressure losses 
        P_loss_viscous = (Loop.f_calc*Loop.rho_c*Loop.vel_c^2*Geo.dl(d))/(2*Loop.D_h_loc);

        if (d ~= length(Geo.pos_i))
            if (Loop.A_c_cs < Loop.A_c_cs_next)
                    K = ((Loop.A_c_cs/Loop.A_c_cs_next)^2-1)^2;
            elseif (Loop.A_c_cs > Loop.A_c_cs_next)
                K = 0.5-0.167*(Loop.A_c_cs_next/Loop.A_c_cs)-...
                    0.125*(Loop.A_c_cs_next/Loop.A_c_cs)^2-...
                    0.208*(Loop.A_c_cs_next/Loop.A_c_cs)^3;
            else
                K = 0;
            end
            P_loss_area = 0.5*K*Loop.rho_c*Loop.vel_c^2;
        else
            P_loss_area = 0;
        end

        P_loss_mom = Param.mdot_f^2*... % Assume den diff is negligible, unless can find a way to get next station den 
            (2/(Loop.A_c_cs*Geo.num_channel+Loop.A_c_cs_next*Geo.num_channel))*...
            (1/(Loop.rho_c*Loop.A_c_cs*Geo.num_channel) - 1/(Loop.rho_c*Loop.A_c_cs_next*Geo.num_channel));

        P_loss_tot = P_loss_mom + P_loss_area + P_loss_viscous;
        Arrays.P_loss_array(d) = P_loss_tot;
        Loop.P_loc = Loop.P_loc - P_loss_tot; 
        Arrays.P_array(d) = convpres(Loop.P_loc, 'Pa', 'psi');
        %Arrays.P_array(d) = Loop.P_loc;

        % Stresses
        Loop.T_iw = (Temp.T_hw + Temp.T_cw)/2;
        Loop.T_ow = Temp.T_ow;
        %Stress = stressAnalysis(Geo, Loop, Mat, Param, d, Temp);
        %Arrays.stress_array(d) = convpres(Stress.sigma_VM, 'Pa', 'psi');
     end


     current_P_error = Loop.P_loc - P_target;
     % Secant
     if Loop.iter_P == 1
         Loop.P_prev_error = current_P_error;
         temp_P = Loop.P_guess;
         Loop.P_guess = Loop.P_prev_guess;
         Loop.P_prev_guess = temp_P;
     else
         P_next = Loop.P_guess - current_P_error * (Loop.P_guess - Loop.P_prev_guess) / (current_P_error - Loop.P_prev_error + 1e-10);
         P_next = max(Pamb * 6894.75, P_next); % Prevent negative or physically impossible pressure guesses
         Loop.P_prev_guess = Loop.P_guess;
         Loop.P_prev_error = current_P_error;
         Loop.P_guess = P_next;
     end
     Loop.P_error = current_P_error;
     if Loop.iter_P > 50
         disp('Pressure loop failed to converge');
         break;
     end
end

%% Plots
fig_configs = {
    'Temperature', {'m', 'r', 'b', 'c', 'g', 'k'}, ...
        {Arrays.T_tc_array, Arrays.T_hw_array, Arrays.T_cw_array, Arrays.T_bulk_array, Arrays.T_fin_array, Arrays.T_sat_array}, ...
        {'Thermal Coating', 'Hot Wall', 'Cold Wall', 'Bulk Coolant', 'Outer Jacket', 'Saturation'}, 'Temperature (K)', 'temperatures.pdf';
    'Pressure', {'b'}, {Arrays.P_array}, {'Coolant Static Pressure'}, 'Pressure (psi)', 'pressure.pdf';
    'HeatFlux', {'r', 'k--', 'b--'}, {Arrays.q_flux_array, Arrays.CHF_array, Arrays.CHF_array_safe}, ...
        {'Heat Flux', 'CHF Limit', 'CHF Limit -10%'}, 'Heat Flux (W/m^2)', 'heatflux.pdf';
    'GasHTC', {'r'}, {Arrays.h_i_array}, {'Gas Heat Transfer Coefficient'}, 'Heat Transfer Coefficient (kg/m^2*s)', 'gashtc.pdf';
    'CoolantHTC', {'b', 'g', 'r', 'm'}, {Arrays.h_c_array, Arrays.h_c_f_array, Arrays.h_nb_array, Arrays.h_tp_array}, ...
        {'Gnielinsky', 'Fin-Corrected Gnielinsky', 'Nucleate Boiling', 'Combined'}, 'Heat Transfer Coefficient (W/m^2*K)', 'coolhtc.pdf'; 
    'Heat Transfer Rate', {'b'}, {Arrays.q_eq_array}, {'Heat Transfer Rate'}, 'Watts', '';
    'Coolant Velocity', {'b'}, {Arrays.vel_c_array}, {'Coolant Velocity'}, 'm/s', '';
};

for i = 1:size(fig_configs, 1)
    figure('Name', fig_configs{i, 1}, 'Color', 'w'); hold on; grid on;
    for j = 1:length(fig_configs{i, 3})
        plot(Geo.pos_i, fig_configs{i, 3}{j}, fig_configs{i, 2}{j}, 'LineWidth', 2, 'DisplayName', fig_configs{i, 4}{j});
    end
    if length(fig_configs{i, 4}) > 1; legend(); end
    title(fig_configs{i, 1}); xlabel('Axial Position x (m)'); ylabel(fig_configs{i, 5});
    xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    if ~isempty(fig_configs{i, 6}); exportgraphics(gcf, fig_configs{i, 6}, 'ContentType','vector'); end
    hold off;
end

%% CSV export
export_list = {
    Geo.pos_i, 'pos_i.csv';
    Gas.pressure, 'P_gas.csv';
    Gas.Taw, 'T_aw.csv';
    Arrays.h_g_array, 'h_g.csv';
    Arrays.P_array, 'P_coolant.csv';
    Arrays.T_bulk_array, 'T_bulk.csv';
    Arrays.h_c_array, 'h_c_fins.csv';
    Arrays.h_c_FEA_array, 'h_c_wall.csv';
    Arrays.gas_flux_array, 'gas_flux.csv';
    Arrays.k_coating_array, 'k_coating.csv';
    };
for y = 1:size(export_list, 1)
    current_data = export_list{y, 1};
    current_filename = export_list{y, 2};

    writematrix(current_data.', current_filename);
end