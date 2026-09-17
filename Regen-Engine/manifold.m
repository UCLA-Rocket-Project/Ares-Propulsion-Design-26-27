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
o_f = 1.25;
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
Geo.num_channel = floor(Geo.circ_t_base / (Geo.w_rib + Geo.min_tol));
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

%% Aft Manifold Geometry
r_channel_base = Geo.D_channel_base ./ 2;
r_outer_jacket = r_channel_base + Geo.h_channel;
r_outer_wall = r_channel_base + Geo.h_channel + Geo.out_wall_thickness;

theta = linspace(0, pi, 100); % Row vector
D_boattail = 0.13025; % m
R_boattail = D_boattail / 2;
R_outer_exit = r_outer_wall(end);
Z_exit = Geo.pos_i(end);

% Downcomer area
A_downcomer = pi * (0.01332)^2; 
A_manifold_inlet = A_downcomer / 2; 
A_manifold_min = Geo.w_channel(end) * Geo.h_channel(end); 

W_max = R_boattail - R_outer_exit;
W_min = Geo.w_channel(end); 

% 1. Target Schedules for Area and Width
A_theta = A_manifold_inlet .* (1 - (theta ./ pi)) + A_manifold_min .* (theta ./ pi);
W_theta = W_max .* (1 - (theta ./ pi)) + W_min .* (theta ./ pi);

% 2. Continuous Height Calculation
% For a conformal wedge, Area approx = 0.5 * Width * Height. 
% To create a bulged outer roof, we use a shape factor (1 = rectangle, 2 = triangle).
shape_factor = 1.6; 
H_raw = shape_factor .* (A_theta ./ W_theta);

% User Correction: Increase terminal height to allow a smooth, manufacturable exit
H_terminal = 0.015; % Set to 15mm to give the blend breathing room

% 3. Apply Tangential Blend to the Height Schedule
blend_percent = 0.75; % Start blending at the 75% mark
blend_idx = floor(blend_percent * length(theta));
H_blend_start = H_raw(blend_idx);

H_smooth = H_raw;
t_blend = linspace(0, 1, length(theta) - blend_idx + 1);
% Quadratic blend targeting the raised terminal height
H_smooth(blend_idx:end) = H_terminal - (H_terminal - H_blend_start) .* (1 - t_blend).^2;

% 4. Map to 3D Coordinates continuously using interp1 (No discrete loops)
Z_top_array = Z_exit - H_smooth;
% Find exact radius on the nozzle wall for every Z-height
R_top_array = interp1(Geo.pos_i, r_outer_wall, Z_top_array, 'linear', 'extrap');

% 5. 3D Cartesian Coordinates for CAD Guide Curves
X_bottom = (R_outer_exit + W_theta) .* cos(theta);
Y_bottom = (R_outer_exit + W_theta) .* sin(theta);
Z_bottom = ones(size(theta)) .* Z_exit;

X_top = R_top_array .* cos(theta);
Y_top = R_top_array .* sin(theta);
Z_top = Z_top_array;

writematrix([X_bottom', Y_bottom', Z_bottom'], 'manifold_guide_bottom.txt');
writematrix([X_top', Y_top', Z_top'], 'manifold_guide_top.txt');

%% Visualization Plot
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

% Visual Bulge Factor (0 = flat line, >0 = convex teardrop roof)
vis_bulge = 0.4;

% Conformal Plotting - Inlet (0 degrees)
Z_roof_0 = linspace(Z_top_array(1), Z_exit, 20);
norm_z_0 = (Z_roof_0 - Z_top_array(1)) ./ (Z_exit - Z_top_array(1));
R_straight_0 = R_top_array(1) + ((R_outer_exit + W_theta(1)) - R_top_array(1)) .* norm_z_0;
R_roof_0 = R_straight_0 + vis_bulge .* ((R_outer_exit + W_theta(1)) - R_top_array(1)) .* norm_z_0 .* (1 - norm_z_0);

idx_0 = find(Geo.pos_i >= Z_top_array(1), 1);
z_inner_0 = [Z_top_array(1), Geo.pos_i(idx_0:end)];
r_inner_0 = [R_top_array(1), r_outer_wall(idx_0:end)];
Z_man_0 = [Z_roof_0, z_inner_0(end:-1:1)];
R_man_0 = [R_roof_0, r_inner_0(end:-1:1)];
plot(Z_man_0, R_man_0, 'm-', 'LineWidth', 2, 'DisplayName', 'Manifold Inlet (0^\circ)');

% Conformal Plotting - End (180 degrees)
Z_roof_180 = linspace(Z_top_array(end), Z_exit, 20);
norm_z_180 = (Z_roof_180 - Z_top_array(end)) ./ (Z_exit - Z_top_array(end));
R_straight_180 = R_top_array(end) + ((R_outer_exit + W_theta(end)) - R_top_array(end)) .* norm_z_180;
R_roof_180 = R_straight_180 + vis_bulge .* ((R_outer_exit + W_theta(end)) - R_top_array(end)) .* norm_z_180 .* (1 - norm_z_180);

idx_180 = find(Geo.pos_i >= Z_top_array(end), 1);
z_inner_180 = [Z_top_array(end), Geo.pos_i(idx_180:end)];
r_inner_180 = [R_top_array(end), r_outer_wall(idx_180:end)];
Z_man_180 = [Z_roof_180, z_inner_180(end:-1:1)];
R_man_180 = [-R_roof_180, -r_inner_180(end:-1:1)]; % Negative for bottom half
plot(Z_man_180, R_man_180, 'm--', 'LineWidth', 2, 'DisplayName', 'Manifold End (180^\circ)');

% Manifold Top Edge Projection
plot(Z_top_array, R_top_array .* cos(theta), 'm-', 'LineWidth', 2, 'DisplayName', 'Manifold Top Edge Profile');
% Exit Plane Vertical Line
plot([Z_exit, Z_exit], [R_outer_exit, -R_outer_exit], 'm', 'LineWidth', 2, 'HandleVisibility', 'off');

title('Regen 1D Profile')
xlabel('Axial Position x (m)');
ylabel('Radial Position y (m)');
axis equal;
xline(0, 'r--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'geometry.pdf', 'ContentType','vector');
hold off;