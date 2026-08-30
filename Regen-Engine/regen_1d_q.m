% Test script for 1d regen engine

clc; clear; close all;
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultTextColor', 'k');

%% Basic Parameters

F = 1450 * 4.44822; % Target thrust in N
Pamb = 12.3; % psia Geo.At 10,500 ft altitude
T_amb = 298; % K (from feed calc sheet)
mfrac_eth = 0.75;
mfrac_h2o = 1 - mfrac_eth;
mw_eth = 46.068; % g/mol
mw_h2o = 18.015; % g/mol
x_eth = (mfrac_eth/mw_eth) / (mfrac_eth/mw_eth + mfrac_h2o/mw_h2o); % Mole fraction
x_h2o = 1 - x_eth;
P_crit = [6140000, 22064000]; % ethanol, water
Param.P_crit_mix = x_eth * P_crit(1) + x_h2o * P_crit(2); % fuel mixture
oxidizer = 'LOX';
% Assumed efficiencies
cstar_eff = 0.90;
cf_eff = 0.98;

% CEA parameters where Pamb = 9.94 psia (interpolate for diff values across nozzle after)
% Taken from throat (A/Geo.At = 1.00)
o_f = 1.3;
Pc_us = 370; % psia, target
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
%Me = sqrt((2/(Gas.gamma-1))*((Pc_us/Pamb)^((Gas.gamma - 1)/Gas.gamma) - 1));
%Geo.Ae = Geo.At * (1/Me)*((2/(Gas.gamma + 1))*(1 + ((Gas.gamma - 1)/2) * Me^2))^((Gas.gamma + 1)/(2*(Gas.gamma - 1)));
Geo.Ae = Geo.At * exp_ratio; % m^2
Geo.R_exit = sqrt(Geo.Ae/pi); % m
Geo.percent_len = 0.8; % input (0.8 is optimal fractional length for most cases)
Geo.len_div = Geo.percent_len * ((Geo.R_exit - Geo.Rt)/tan(deg2rad(15))); % m diverging length from radii
Geo.theta_e = deg2rad(13); % deg, from HH fig 4.16
Geo.theta_n = deg2rad(23); % deg

%% Chamber and Converging Geometry
Geo.conv_angle = deg2rad(45); % deg, standard
Geo.L_star = 40; % in, optimal for ethanol/lox
Geo.CR = 5; % A_c / A_t
%Geo.id_chanber = 4.75; % in, heritage
Geo.V_total = convlength(Geo.L_star, 'in', 'm') * Geo.At; % m^3
%Geo.Rc = convlength(Geo.id_chanber / 2, 'in', 'm'); % m, chamber radius
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

% Recentering coordinates to around throat
Geo.x_t = 0;
Geo.x_chamber_start = -(Geo.len_chamber + Geo.len_conv);
Geo.x_exit = Geo.len_div;
% Position arrrays
Geo.pos_i = Geo.x_chamber_start:Geo.dx:Geo.x_exit; % axial position
Geo.pos_j = zeros(size(Geo.pos_i)); % radii
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
for k = 1:length(Geo.pos_i)
    x = Geo.pos_i(k);
    if x < Geo.x_chamber_end % chamber
        Geo.pos_j(k) = Geo.Rc;
    elseif x >= Geo.x_chamber_end && x < Geo.x_f_tangent
        Geo.pos_j(k) = Geo.y_conv_fillet + sqrt(Geo.R_fillet^2 - (x - Geo.x_conv_fillet)^2);
    elseif x >= Geo.x_f_tangent && x < Geo.x_conv_tangent % straight converging cone
        Geo.pos_j(k) = Geo.y_f_tangent - tan(Geo.conv_angle) * (x - Geo.x_f_tangent);
    elseif x >= Geo.x_conv_tangent && x < 0 % converging throat arc
        Geo.pos_j(k) = (Geo.Rt + Geo.R_curve) - sqrt((Geo.R_curve)^2 - x^2); % Geo.Rt + arc length - height of curve Geo.At point
    elseif x >= 0 && x < Geo.xn % diverging throat arc
        Geo.pos_j(k) = (Geo.Rt + 0.382*Geo.Rt) - sqrt((0.382*Geo.Rt)^2 - x^2);
    else % Rao parabola approx
        Geo.pos_j(k) = Geo.a_Rao * x^2 + b_Rao * x + c_Rao;
    end
end
% Slant Length
Geo.dx_slope = gradient(Geo.pos_j, Geo.dx);
Geo.dl = Geo.dx .* sqrt(1 + Geo.dx_slope .^2);

%% Cooling Channel and Outer Jacket Geometry
% Fixed channel height and wall thickness
Geo.min_tol = 0.001; % m 3d printer tolerance
Geo.D_gas = 2 * Geo.pos_j; % Array of gas-side diameter Geo.At every node
Geo.D_t = 2 * Geo.Rt;
Geo.h_channel = Geo.min_tol*2; % channel height (radial)
Geo.wall_thickness = Geo.min_tol; % HW/Loop.cw
Geo.w_rib = Geo.min_tol; % fixed, rib width
Geo.D_channel_base = Geo.D_gas + 2 * Geo.wall_thickness; % Engine diameters with added wall thickness
Geo.D_t_base = Geo.D_t + 2 * Geo.wall_thickness; % Throat diameter with added wall thickness
% Number of channels determined Geo.At throat
Geo.circ_t_base = pi * (Geo.D_t_base); % Gas side diameter + wall thickness
Geo.num_channel = floor(Geo.circ_t_base / (Geo.w_rib + Geo.min_tol));
% Local circumferences across engine
Geo.circ_local_base = pi * Geo.D_channel_base;
Geo.w_channel = (Geo.circ_local_base - (Geo.num_channel * Geo.w_rib)) ./ Geo.num_channel; % variable, channel widths
Geo.D_h = (4 .* Geo.w_channel .* Geo.h_channel) ./ (2 * Geo.w_channel + 2 * Geo.h_channel); % hydraulic diameter of rectangular channels
Geo.Per_heated = Geo.w_channel + 2 * Geo.h_channel; % heated perimeter

%% HT Areas
Geo.A_gas = pi .* Geo.D_gas .* Geo.dl; % Gas-wall convection SA
Geo.A_w = pi .* ((Geo.D_gas + Geo.D_channel_base)./2) .* Geo.dl; % wall-wall conduction SA (average diameter)
Geo.A_co = Geo.Per_heated .* Geo.num_channel .* Geo.dl; % coolant side surface area (heated)
Geo.A_wc = Geo.w_channel .* Geo.dl; % cool wall area per increment

%% Visualization Plot
figure('Name', '1D Engine Geometry', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Geo.pos_j, 'k', 'LineWidth', 2, 'DisplayName', 'Hot Wall');
plot(Geo.pos_i, -Geo.pos_j, 'k', 'LineWidth', 2, 'HandleVisibility','off');

r_channel_base = Geo.D_channel_base ./ 2;
plot(Geo.pos_i, r_channel_base, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Cold Wall');
plot(Geo.pos_i, -r_channel_base, 'b--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

r_outer_jacket = r_channel_base + Geo.h_channel;
plot(Geo.pos_i, r_outer_jacket, 'b', 'LineWidth', 2, 'DisplayName', 'Outer Jacket');
plot(Geo.pos_i, -r_outer_jacket, 'b', 'LineWidth', 2, 'HandleVisibility', 'off');

title('Regen 1D Profile')
xlabel('Axial Position x (m)');
ylabel('Radial Position y (m)');
axis equal;
xline(0, 'r--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'geometry.pdf', 'ContentType','vector');
hold off;


%% Gas Properties
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

Gas.pressure = Param.Pc * (1 + (Gas.gamma - 1)/2 .* Gas.M_local .^2) .^ (-Gas.gamma / (Gas.gamma - 1));
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

%% Coolant Properties
% Required properties update for equilibrium: Cp, k, mu, rho
% Define Table bounds (Tmax < Loop.T_sat Geo.At Pmin or coolprop will crash)
table_filename = 'coolprop_tables.mat';
table_version_req = 3; % Bump to force cached tables to rebuild
tables_loaded = false;
if isfile(table_filename)
    tbl = load(table_filename);
    if isfield(tbl, 'table_version') && tbl.table_version == table_version_req
        %get_rho = tbl.get_rho; get_cp = tbl.get_cp; get_mu = tbl.get_mu; get_k = tbl.get_k;
        %get_T_sat = tbl.get_T_sat; get_rho_l = tbl.get_rho_l; get_rho_v = tbl.get_rho_v;
        %get_surften = tbl.get_surften; get_h_fg = tbl.get_h_fg; get_P_sat = tbl.get_P_sat;
        P_vec = tbl.P_vec; T_vec = tbl.T_vec;
        Cool = tbl.Cool;
        tables_loaded = true;
    end
end
if ~tables_loaded % create tables first time (rebuild when table_version_req is bumped)
    res = 50; % 50x50 data grid
    P_min = convpres(300, 'psi', 'Pa'); % Pa, Geo.At chamber
    P_max = convpres(900, 'psi', 'Pa'); % Pa, Geo.At fuel tank
    T_min = 290; % K, standard inlet temp
    T_max = 450; % K, below 300 psi boiling point of ethanol
    % 1D vectors for P and T (fast solve for Coolprop, easy to get)
    P_vec = linspace(P_min, P_max, res);
    T_vec = linspace(T_min, T_max, res);

    A12_vl = 1.6798; A21_vl = 0.9227; % Van Laar coeffs
    % Activity coefficients from Van Laar equation: deviation of a mixture of chemical substances from ideal behaviour
    gam_eth = exp(A12_vl*(A21_vl*x_h2o/(A12_vl*x_eth + A21_vl*x_h2o))^2); 
    gam_h2o = exp(A21_vl*(A12_vl*x_eth/(A12_vl*x_eth + A21_vl*x_h2o))^2);

    % use ndgrid for higher dimensionality use
    [P_grid, T_grid] = ndgrid(P_vec, T_vec);

    % 2D grids for bulk properties
    rho_grid = zeros(res,res);
    cp_grid = zeros(res,res);
    mu_grid = zeros(res,res);
    k_grid = zeros(res,res);
    % 1D grids
    T_sat_grid = zeros(1,res);
    rho_l_grid = zeros(1,res); % sat. liquid density
    rho_v_grid = zeros(1,res); % sat. vapor density
    surften_grid = zeros(1,res); % surface tension
    h_fg_grid = zeros(1,res); % latent heat of vaporization

    % Generate values
    for i = 1:res
        P_val = P_vec(i);
        T_sat_h2o = py.CoolProp.CoolProp.PropsSI('T','P',P_val,'Q',0, 'water');
        if P_val < 6100000 % Ethanol critical pressure in Pa (with buffer)
            T_sat_eth = py.CoolProp.CoolProp.PropsSI('T','P',P_val,'Q',0, 'ethanol');
            % T_sat bisec iteration based on Raoult's Law to find bubble point
            T_sat_lo = T_sat_eth - 15; 
            T_sat_hi = min(T_sat_h2o, 513.5); % Capped at ethanol critical temperature
            for t = 1:40
                T_sat_mid = 0.5*(T_sat_lo + T_sat_hi);
                P_bubble = x_eth*gam_eth*py.CoolProp.CoolProp.PropsSI('P','T',T_sat_mid,'Q',0,'ethanol') + ...
                    x_h2o*gam_h2o*py.CoolProp.CoolProp.PropsSI('P','T',T_sat_mid,'Q',0,'water');
                if P_bubble < P_val
                    T_sat_lo = T_sat_mid;
                else
                    T_sat_hi = T_sat_mid;
                end
            end
            T_sat_grid(i) = 0.5*(T_sat_lo + T_sat_hi);
            % Liquid Density
            rho_l_eth = py.CoolProp.CoolProp.PropsSI('D','P',P_val,'Q',0, 'ethanol');
            rho_l_h2o = py.CoolProp.CoolProp.PropsSI('D','P',P_val,'Q',0, 'water');
            rho_l_grid(i) = 1 / ((mfrac_eth / rho_l_eth) + (mfrac_h2o / rho_l_h2o));
            % Vapor Density
            rho_v_eth = py.CoolProp.CoolProp.PropsSI('D','P',P_val,'Q',1, 'ethanol');
            rho_v_h2o = py.CoolProp.CoolProp.PropsSI('D','P',P_val,'Q',1, 'water');
            rho_v_grid(i) = 1 / ((mfrac_eth / rho_v_eth) + (mfrac_h2o / rho_v_h2o));
            % Surface Tension (Interfacial)
            surften_eth = py.CoolProp.CoolProp.PropsSI('I','P',P_val,'Q',0, 'ethanol');
            surften_h2o = py.CoolProp.CoolProp.PropsSI('I','P',P_val,'Q',0, 'water');
            surften_grid(i) = mfrac_eth * surften_eth + mfrac_h2o * surften_h2o;
            % Latent Heat of Vaporization (H_vap - H_liq)
            h_v_eth = py.CoolProp.CoolProp.PropsSI('H','P',P_val,'Q',1, 'ethanol');
            h_l_eth = py.CoolProp.CoolProp.PropsSI('H','P',P_val,'Q',0, 'ethanol');
            h_fg_eth = h_v_eth - h_l_eth;
            h_v_h2o = py.CoolProp.CoolProp.PropsSI('H','P',P_val,'Q',1, 'water');
            h_l_h2o = py.CoolProp.CoolProp.PropsSI('H','P',P_val,'Q',0, 'water');
            h_fg_h2o = h_v_h2o - h_l_h2o;
            h_fg_grid(i) = mfrac_eth * h_fg_eth + mfrac_h2o * h_fg_h2o;
        else % undefined for superheated vapor
            T_sat_grid(i) = NaN;
            rho_l_grid(i) = NaN;
            rho_v_grid(i) = NaN;
            surften_grid(i) = NaN;
            h_fg_grid(i) = NaN;
        end
        for j = 1:res
            P_val = P_grid(i,j);
            T_val = T_grid(i,j);
            
            rho_eth = py.CoolProp.CoolProp.PropsSI('D','T',T_val,'P',P_val,'ethanol');
            rho_h2o = py.CoolProp.CoolProp.PropsSI('D','T',T_val,'P',P_val,'water');
            cp_eth = py.CoolProp.CoolProp.PropsSI('C','T',T_val,'P',P_val,'ethanol');
            cp_h2o = py.CoolProp.CoolProp.PropsSI('C','T',T_val,'P',P_val,'water');
            mu_eth = py.CoolProp.CoolProp.PropsSI('V','T',T_val,'P',P_val,'ethanol');
            mu_h2o = py.CoolProp.CoolProp.PropsSI('V','T',T_val,'P',P_val,'water');
            k_eth = py.CoolProp.CoolProp.PropsSI('L','T',T_val,'P',P_val,'ethanol');
            k_h2o = py.CoolProp.CoolProp.PropsSI('L','T',T_val,'P',P_val,'water');
            
            rho_grid(i,j) =  1 / ((mfrac_eth / rho_eth) + (mfrac_h2o / rho_h2o));
            cp_grid(i,j) = mfrac_eth * cp_eth + mfrac_h2o * cp_h2o;
            % Grunberg-Nissan relation for viscosity of fluid mixture
            G12 = 840 / T_val;
            mu_grid(i,j) = exp(x_eth*log(mu_eth) + x_h2o*log(mu_h2o) + x_eth*x_h2o*G12);
            k_grid(i,j) = mfrac_eth * k_eth + mfrac_h2o * k_h2o;
        end
    end
    % Interpolation objects
    Cool.get_rho = griddedInterpolant(P_grid, T_grid, rho_grid, 'linear', 'nearest');
    Cool.get_cp = griddedInterpolant(P_grid, T_grid, cp_grid, 'linear', 'nearest');
    Cool.get_mu = griddedInterpolant(P_grid, T_grid, mu_grid, 'linear', 'nearest');
    Cool.get_k = griddedInterpolant(P_grid, T_grid, k_grid, 'linear', 'nearest');
    sat_idx = ~isnan(T_sat_grid);
    Cool.get_T_sat = griddedInterpolant(P_vec(sat_idx), T_sat_grid(sat_idx), 'linear', 'nearest');
    Cool.get_rho_l = griddedInterpolant(P_vec(sat_idx), rho_l_grid(sat_idx), 'linear', 'nearest');
    Cool.get_rho_v = griddedInterpolant(P_vec(sat_idx), rho_v_grid(sat_idx), 'linear', 'nearest');
    Cool.get_surften = griddedInterpolant(P_vec(sat_idx), surften_grid(sat_idx), 'linear', 'nearest');
    Cool.get_h_fg = griddedInterpolant(P_vec(sat_idx), h_fg_grid(sat_idx), 'linear', 'nearest');
    Cool.get_P_sat = griddedInterpolant(T_sat_grid(sat_idx), P_vec(sat_idx), 'linear', 'nearest');

    table_version = table_version_req;

    save(table_filename,'Cool','P_vec', 'T_vec', 'table_version');
end

%% Material Properties (316 Stainless)
T_mp = 1643.15; % K
% Wall thermal conductivity k_w = f(T)
Mat.k_w_ref_temps = [-0.15, 19.85, 26.85, 76.85, 126.85, 226.85, 326.85, 426.85, 526.85, 626.85, ...
    726.85, 826.85, 926.85, 1026.9, 1126.9, 1226.9, 1326.9, 1370.9, 1398.9, 1426.9] + 273.15; % K
Mat.k_w_ref = [12.97, 13.31, 13.44, 14.32, 15.16, 16.8, 18.36, 19.87, 21.39, 22.79, 24.06, 25.46, 26.74, ...
    28.02, 29.32, 30.61, 31.86, 32.41, 26.9, 27.24];
Mat.r = 1e-4; % m, surface roughness
%get_k_w = interp1(Mat.k_w_ref_temps, Mat.k_w_ref, T_w_local, 'linear'); % callout

%% Output Arrays
Arrays.T_hw_array = zeros(size(Geo.pos_i));
Arrays.T_cw_array = zeros(size(Geo.pos_i));
Arrays.T_bulk_array = zeros(size(Geo.pos_i));
Arrays.P_array = zeros(size(Geo.pos_i));
Arrays.P_loss_array = zeros(size(Geo.pos_i));
Arrays.q_flux_array = zeros(size(Geo.pos_i));
Arrays.h_i_array = zeros(size(Geo.pos_i));
Arrays.h_c_array = zeros(size(Geo.pos_i));
Arrays.h_c_f_array = zeros(size(Geo.pos_i));
Arrays.h_nb_array = zeros(size(Geo.pos_i));
Arrays.CHF_array = zeros(size(Geo.pos_i));
Arrays.stress_array = zeros(size(Geo.pos_i));
Arrays.q_eq_array = zeros(size(Geo.pos_i));

% Test Arrays (remove when done)
Arrays.f_array = zeros(size(Geo.pos_i));
Arrays.Nu_array = zeros(size(Geo.pos_i));
Arrays.T_sat_array = zeros(size(Geo.pos_i));
Arrays.Re_array = zeros(size(Geo.pos_i));
Arrays.vel_c_array = zeros(size(Geo.pos_i));
Arrays.fin_eff_array = zeros(size(Geo.pos_i));
Arrays.q_error_array = zeros(size(Geo.pos_i));
Arrays.iter_T_array = zeros(size(Geo.pos_i));
Arrays.sigma_array = zeros(size(Geo.pos_i));

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
    Loop.T_bulk = T_amb;

    for d = length(Geo.pos_i):-1:1 % Axial marching loop
        % Local Geometry, add channel height array earlier
        Loop.cw = Geo.w_channel(d);
        Loop.ch = Geo.h_channel;
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
        Arrays.T_cw_array(d) = Temp.T_cw;
  
        Arrays.T_hw_array(d) = Temp.T_hw_guess;
        Arrays.fin_eff_array(d) = Temp.fin_eff;
        Arrays.h_i_array(d) = Temp.h_i;
        Arrays.q_error_array(d) = Temp.q_error;
        Arrays.iter_T_array(d) = Temp.iter_T;
        Arrays.q_eq_array(d) = Temp.q_eq;
        Arrays.h_c_f_array(d) = Temp.h_c_f;
        Arrays.h_nb_array(d) = Temp.h_nb;
        Arrays.sigma_array(d) = Temp.sigma;
        
        % Check CHF
        q_flux = Temp.q_eq/Loop.A_conv;
        Arrays.q_flux_array(d) = q_flux;
        dT_sub = Loop.T_sat - Loop.T_bulk; % K, bulk subcooling
        F_p = 1.17-8.56*(10^(-4))*convpres(Loop.P_loc, 'Pa', 'psi');
        CHF_base = 0.1003+0.05264*sqrt(convvel(Loop.vel_c, 'm/s', 'ft/s')*(dT_sub * 9/5)); 
        Arrays.CHF_array(d) = CHF_base*F_p*1635000;
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

        % Calculate Von Mises Stress
       
        % i disabled for single loop, needs Loop.T_close_ow/T_far_ow, 
        % Geo.out_wall_thickness, Mat.alpha/E/nu defined (Loop.T_cw was renamed Temp.T_cw) ----- manav
        
        %Loop.T_iw = (Temp.T_hw_guess - Temp.T_cw)/2;
        %Loop.T_ow = (Loop.T_close_ow - Loop.T_far_ow);
        %stressAnalysis(Geo, Loop, Mat, Param, d);
        %Arrays.stress_array(d) = Stress.sigma_VM;
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
         % Prevent negative or physically impossible pressure guesses
         P_next = max(Pamb * 6894.75, P_next);
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

% Temps
figure('Name', 'Temperature', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.T_hw_array, 'r', 'LineWidth', 2);
plot(Geo.pos_i, Arrays.T_cw_array, 'b', 'LineWidth', 2);
plot(Geo.pos_i, Arrays.T_bulk_array, 'c', 'LineWidth', 2);
plot(Geo.pos_i, Arrays.T_sat_array, 'g', 'LineWidth', 2);
legend('Hot Wall', 'Cold Wall', 'Bulk Coolant')
title('Temperatures')
xlabel('Axial Position x (m)');
ylabel('Temperature (K)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'temperatures.pdf', 'ContentType','vector');
hold off;

% Pressure
figure('Name', 'Pressure', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.P_array, 'b', 'LineWidth', 2);
title('Coolant Static Pressure')
xlabel('Axial Position x (m)');
ylabel('Pressure (psi)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'pressure.pdf', 'ContentType','vector');

hold off;

% Heat Flux
figure('Name', 'HeatFlux', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.q_flux_array, 'r', 'LineWidth', 2);
plot(Geo.pos_i, Arrays.CHF_array, 'k--', 'LineWidth', 1.5, 'DisplayName', 'CHF Limit');
title('Coolant Heat Flux')
xlabel('Axial Position x (m)');
ylabel('Heat Flux (W/m^2)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'heatflux.pdf', 'ContentType','vector');

hold off;

% HTC
figure('Name', 'GasHTC', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.h_i_array, 'r', 'LineWidth', 2);
title('Gas Heat Transfer Coefficient')
xlabel('Axial Position x (m)');
ylabel('Heat Transfer Coefficient (W/m^2*K)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'gashtc.pdf', 'ContentType','vector');
hold off;

figure('Name', 'CoolantHTC', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.h_c_array, 'b', 'LineWidth', 2);
plot(Geo.pos_i, Arrays.h_c_f_array, 'g', 'LineWidth', 2);
plot(Geo.pos_i, Arrays.h_nb_array, 'r', 'LineWidth', 2);
legend('Gnielinsky', 'Fin-Corrected Gnielinsky', 'Nucleate Boiling')
title('Coolant Heat Transfer Coefficient')
xlabel('Axial Position x (m)');
ylabel('Heat Transfer Coefficient (W/m^2*K)');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'coolhtc.pdf', 'ContentType','vector');
hold off;

% Q
figure('Name', 'Heat Transfer Rate', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.q_eq_array, 'b', 'LineWidth', 2);
title('Heat Transfer Rate')
xlabel('Axial Position x (m)');
ylabel('Watts');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

%{
% Friction Factor
figure('Name', 'Darcy Friction Factor', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.f_array, 'b', 'LineWidth', 2);
title('Friction Factor')
xlabel('Axial Position x (m)');
ylabel('Dimensionless');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Nu
figure('Name', 'Nusselt Number', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.Nu_array, 'b', 'LineWidth', 2);
title('Nusselt Number')
xlabel('Axial Position x (m)');
ylabel('Dimensionless');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Saturation Temp
figure('Name', 'Saturation Temp', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.T_sat_array, 'b', 'LineWidth', 2);
title('Saturation Temp')
xlabel('Axial Position x (m)');
ylabel('K');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Re
figure('Name', 'Geo.Re', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.Re_array, 'b', 'LineWidth', 2);
title('Geo.Re')
xlabel('Axial Position x (m)');
ylabel('Dimensionless');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Coolant Velocity
figure('Name', 'Coolant Velocity', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.vel_c_array, 'b', 'LineWidth', 2);
title('Coolant Velocity')
xlabel('Axial Position x (m)');
ylabel('m/s');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Fin Eff
figure('Name', 'Fin Efficiency', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.fin_eff_array, 'b', 'LineWidth', 2);
title('Fin Efficiency')
xlabel('Axial Position x (m)');
ylabel('Dimensionless');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Q Error
figure('Name', 'Q error', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.q_error_array, 'b', 'LineWidth', 2);
title('Q Error')
xlabel('Axial Position x (m)');
ylabel('Watts');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Iterations per Station
figure('Name', 'Iterations per Station', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.iter_T_array, 'b', 'LineWidth', 2);
title('Iterations per Station')
xlabel('Axial Position x (m)');
ylabel('Number of Iterations');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Sigma
figure('Name', 'Sigma', 'Color', 'w');
hold on; grid on;
plot(Geo.pos_i, Arrays.sigma_array, 'k', 'LineWidth', 2);
title('Sigma')
xlabel('Axial Position x (m)');
ylabel('Dimensionless');
xline(0, 'k--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

%}


%% Thermal Logic

function Temp = temp_iteration(Param, Cantera, Y_str, Geo, Gas, Cool, Mat, Loop, d) % HW Temp Iteration
    Temp.T_hw_guess = 700; % Initial Guess (1.25 FOS applied to material melting point)
    T_hw_prev_guess = 650; % Cooler lower bound (Luca's value)
    tol_q = 0.1; % w
    q_prev_error = 0;
    Temp.iter_T = 0;
    while true
        Temp.iter_T = Temp.iter_T + 1;
        % Fin Efficiency
        Temp.k_w_loc = interp1(Mat.k_w_ref_temps, Mat.k_w_ref, Temp.T_hw_guess, 'linear', 'extrap');
        fin_m = sqrt((2*Loop.h_c)/(Temp.k_w_loc * Geo.w_rib));
        Temp.fin_eff = tanh(fin_m * Loop.ch) / (fin_m * Loop.ch);
        Temp.h_c_f = Loop.h_c*(Loop.cw+2*Temp.fin_eff*Loop.ch)/(Loop.cw+Geo.w_rib); % Fin corrected Loop.h_c

        % Gas convection HTC with Bartz
        Temp.sigma = 1 / ...
            ((0.5 * (Temp.T_hw_guess/Gas.T_stag) * (1 + (Gas.gamma-1)/2 * Gas.M_local(d)^2) + 0.5)^0.68 *...
            (1 + (Gas.gamma-1)/2 * Gas.M_local(d)^2)^0.12);
        
        Temp.h_g = ((0.026/Geo.D_t^0.2)*...
            ((Gas.mu_g_local(d)^0.2*Gas.cp_g_local(d))/Gas.prandtl_g_local(d)^0.6)*...
            (Param.Pc/Param.cstar_act)^0.8)*...
            (Geo.D_t/Geo.R_curve)^0.1*...
            (Geo.At/Gas.A_local(d))^0.9;

        Temp.h_i = ((0.026/Geo.D_t^0.2)*...
            (Gas.mu_g_local(d)^0.2/Gas.prandtl_g_local(d)^0.6)*...
            (Param.Pc/Param.cstar_act)^0.8)*...
            (Geo.D_t/Geo.R_curve)^0.1*...
            (Geo.At/Gas.A_local(d))^0.9*...
            Temp.sigma;

        %Temp.q_eq = Temp.h_g*Loop.A_g_loc*(Loop.Taw_loc - Temp.T_hw_guess);
        
        i_w = get_iw(Cantera, Y_str, Gas, Temp);
        Temp.q_eq = Temp.h_i * Loop.A_g_loc * (Gas.i_aw(d) - i_w);
        Temp.T_cw = Temp.T_hw_guess - (Temp.q_eq)*...
            log(1+2*Geo.wall_thickness/Loop.D_g_loc)/(2*pi*Geo.dl(d)*Temp.k_w_loc); % Cold wall temp derived from guess

        % coolant side shares the convection term in both regimes so q_c is continuous at t_cw = t_sat
        mass_flux = Param.mdot_f / Loop.A_conv;
        dT_sub = Loop.T_sat - Loop.T_bulk;
        dT_bulk = Temp.T_cw - Loop.T_bulk;
        dT_sat = Temp.T_cw - Loop.T_sat;
        Stanton = abs(Temp.q_eq/Loop.A_base_loc) / (mass_flux * Loop.cp_c * dT_sub); % dimensionless
        Peclet = mass_flux * Loop.D_h_loc * Loop.cp_c / Loop.k_c;
        Boiling = abs(Temp.q_eq/Loop.A_base_loc) / (mass_flux * Loop.h_fg);

        % Elfaham Correlation
        
        % Find Ms correction term
        if 1e-10 <= Boiling && Boiling < 1e-3 % actual lower bound 1e-5, change to this if Boiling within appropriate range
            Ms = 0.7;
        elseif 1e-3 <= Boiling && Boiling < 5e-3
            Ms = 1.5;
        elseif 5e-3 <= Boiling && Boiling < 1e-2
            Ms = 1.3;
        else
            Ms = 1.1;
        end
        q_c = Temp.h_c_f*Loop.A_base_loc*(Temp.T_cw - Loop.T_bulk);
        Temp.h_nb = 0;
        if Temp.T_cw > Loop.T_sat
            Temp.h_nb = 55 * (abs(Temp.q_eq/Loop.A_base_loc))^0.67 * Loop.P_reduced^0.12 * (-log10(Loop.P_reduced))^(-0.55) * Param.MW^(-0.5);
            E = (1 + Loop.quality*Loop.prandtl_c*(Loop.rho_c_l/Loop.rho_c_v - 1))^0.35;
            S = Ms / (1 + 0.055 * E^0.1 * Loop.Re^0.16);
            Temp.h_tp = sqrt((S*Temp.h_nb)^2 + (E*Temp.h_c_f)^2);
            q_c = Temp.h_tp * Loop.A_base_loc*(Temp.T_cw - Loop.T_bulk);
        end

        % Zhu-Bi-Yan Correlation        
        % if Stanton <= 38*Peclet^(-0.38) % partially boiling
        %     Temp.h_nb = 0.00122*...
        %         (((Loop.k_c^0.79)*(Loop.cp_c^0.45)*(Loop.rho_c_l^0.49))/...
        %         ((Loop.surften^0.5)*(Loop.mu_c^0.29)*(Loop.h_fg^0.24)*(Loop.rho_c_v^0.24)))*...
        %         ((dT_sat)^0.24)*...
        %         (max(0, Cool.get_P_sat(min(Temp.T_cw, 513)) - Loop.P_loc))^0.75; % Cap at critical pressure for ethanol
        %     S = (1 / (1 + 0.055*Loop.Re^0.16)) * (Loop.T_sat/dT_sub)^0.28;
        %     Temp.h_tp = sqrt(Temp.h_c_f^2 + (S*Temp.h_nb*dT_sat/dT_bulk)^2);
        %     q_c = Temp.h_tp * Loop.A_base_loc*(Temp.T_cw - Loop.T_bulk);
        % else % fully developed boiling
        %     q_c = 1000 * Loop.A_base_loc * dT_sat / (32 * e^(-Loop.P_loc/8.6*10^6)); % W, Zhu-Bi
        % end
        

        % Chen Correlation
        %{
        q_c = Temp.h_c_f*Loop.A_base_loc*(Temp.T_cw - Loop.T_bulk);
        if Temp.T_cw > Loop.T_sat
            Temp.h_nb = 0.00122*...
                (((Loop.k_c^0.79)*(Loop.cp_c^0.45)*(Loop.rho_c_l^0.49))/...
                ((Loop.surften^0.5)*(Loop.mu_c^0.29)*(Loop.h_fg^0.24)*(Loop.rho_c_v^0.24)))*...
                ((Temp.T_cw-Loop.T_sat)^0.24)*...
                (max(0, Cool.get_P_sat(min(Temp.T_cw, 513)) - Loop.P_loc))^0.75; % Cap at critical pressure for ethanol
            S = 1/(1+(2.53*(10^-6))*(Loop.Re^1.17));
            q_c = q_c + S*Temp.h_nb*(Temp.h_c_f/Loop.h_c)*Loop.A_base_loc*(Temp.T_cw - Loop.T_sat);
        end
        %}

        % Secant
        Temp.q_error = Temp.q_eq - q_c;
        if abs(Temp.q_error) < tol_q
            break;
        end
        if Temp.iter_T > 100
            warning('t_hw failed to converge at station %d (residual %.3g w)', d, Temp.q_error);
            break;
        end
        if Temp.iter_T == 1
            q_prev_error = Temp.q_error;
            temp_T = Temp.T_hw_guess;
            Temp.T_hw_guess = T_hw_prev_guess;
            T_hw_prev_guess = temp_T;
        else
            T_next = Temp.T_hw_guess - Temp.q_error * (Temp.T_hw_guess - T_hw_prev_guess) / (Temp.q_error - q_prev_error + 1e-10); % small buffer to avoid divide by 0
            T_next = max(Loop.T_bulk + 1, min(T_next, Loop.Taw_loc - 1));
            T_hw_prev_guess = Temp.T_hw_guess;
            q_prev_error = Temp.q_error;
            Temp.T_hw_guess = T_next;
        end
    end

    function i_w = get_iw(Cantera,Y_str, Gas, Temp)

    Cantera.TPY = py.tuple({Temp.T_hw_guess, Gas.pressure(d), Y_str});
    Cantera.equilibrate('TP');
    i_w = double(Cantera.enthalpy_mass);

    end

end

%% Stresses

function stressAnalysis(Geo, Loop, Mat, Param, d)
    p_diff = Param.Pc - Loop.P_loc;
    r_in = Loop.D_g_loc/2;
    r_out = r_in+Geo.wall_thickness+Loop.ch+Geo.out_wall_thickness;
    r_avg = (r_in + r_out)/2;

    % Shear
    Stress.sigma_s_mech = 0.5*p_diff*Loop.cw*Geo.wall_thickness;
    Stress.sigma_s_therm = ((Loop.cw+Geo.w_rib)*Mat.alpha*Mat.E*...
        (Loop.T_iw-Loop.T_ow)*Geo.wall_thickness*Geo.out_wall_thickness)/...
        (5*Geo.w_rib*(1-Mat.nu)*...
        (Geo.wall_thickness+Geo.out_wall_thickness)^2);
    Stress.sigma_s_tot = Stress.sigma_s_mech + Stress.sigma_s_therm; %need E, nu, temps, out thick

    % Hoop
    Stress.sigma_h_mech = (p_diff * r_avg)/...
        (Geo.wall_thickness+Geo.out_wall_thickness);
    Stress.sigma_h_therm = (Mat.alpha*Mat.E*(Loop.T_iw-Loop.T_ow)*...
        Geo.out_wall_thickness) / ((1-Mat.nu)*...
        (Geo.wall_thickness + Geo.out_wall_thickness));
    Stress.sigma_h_tot = Stress.sigma_h_mech + Stress.sigma_h_therm;

    % Bending
    Stress.sigma_b = (p_diff * Loop.cw^2)/(2*Geo.wall_thickness);
    
    % Axial
    if isnan(d) %less than whatever the throat station is
        Stress.sigma_a = (Gas.P_exit*pi*((Loop.D_g_loc/2)^2-Loop.Geo.Rt^2))/...
            ((pi*(r_out^2 - Loop.Geo.Rt^2))-...
            (Geo.num_channel*Loop.cw*Loop.ch));
    else
        Stress.sigma_a = 0;
    end

    Stress.sigma_VM = sqrt(0.5*((Stress.sigma_a-Stress.sigma_h_tot)^2+...
        (Stress.sigma_h_tot - Stress.sigma_b)^2 +...
        (Stress.sigma_b - Stress.sigma_a)^2) + 3*Stress.sigma_s_tot^2);
end



