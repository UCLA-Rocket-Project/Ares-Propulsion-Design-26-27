clear; clc; close all;

%% INPUTS %%

margin_psi = 100; % psi

P_ambient_psi = 14.7; % psi
P_prefire_psi = 535; % psi
P_tank_burn_psi = 488; % psi

OF = 1.30;
prop_mass_lbm = 100; % lbm
burn_time = 15.06; % s

V_lox_tank_L = 25.0; % L
V_fuel_tank_L = 25.0; % L
V_copv_L = 12.0; % L

rho_lox = 1129.6; % kg/m^3
rho_fuel = 853.4; % kg/m^3 (75 wt% ethanol, 25 wt% water at 480 psi, which is burn tank pressure)

P_copv_start_psi = 4500; % pressure BEFORE ISO VALVE OPENS (tanks at ambien), psi
T_copv_start = 293.15; % K
polytropic_n = 1.667; % COPV polytropic exponent, assume isentropic process, for He, gamma=1.667

T_ullage_initial = 293.15; % K
T_ullage_prefire = 293.15; % sketchy, K

dt = 0.01; % timestep, s


%% CONVERSIONS %%

PSI_TO_PA = 6894.757293168;
LBM_TO_KG = 0.45359237;
L_TO_M3 = 1e-3;
M3S_TO_CFM = 2118.880003;

R_He = 2077.264; % helium gas constant, J/(kg K)
gamma_He = 1.667;
cv_He = R_He / (gamma_He - 1); % specific heat at constant volume
cp_He = gamma_He * cv_He; % specific heat at constant pressure

P_ambient = P_ambient_psi * PSI_TO_PA;
P_prefire = P_prefire_psi * PSI_TO_PA;
P_tank_burn = P_tank_burn_psi * PSI_TO_PA;
P_copv = P_copv_start_psi * PSI_TO_PA;

V_copv = V_copv_L * L_TO_M3;
V_lox_tank = V_lox_tank_L * L_TO_M3;
V_fuel_tank = V_fuel_tank_L * L_TO_M3;

prop_mass = prop_mass_lbm * LBM_TO_KG;

%% INITIAL CALCULATIONS %%
m_fuel = prop_mass / (1 + OF);
m_lox = prop_mass - m_fuel;

mdot_total = prop_mass / burn_time;
mdot_fuel = mdot_total / (1 + OF);
mdot_lox = mdot_total - mdot_fuel;

V_lox_liquid = m_lox / rho_lox;
V_fuel_liquid = m_fuel / rho_fuel;

V_ullage_lox = V_lox_tank - V_lox_liquid;
V_ullage_fuel = V_fuel_tank - V_fuel_liquid;

if V_ullage_lox <= 0
    error('LOX volume is bigger than LOX tank volume.');
end

if V_ullage_fuel <= 0
    error('Fuel volume is bigger than fuel tank volume.');
end

%% INITIAL STATE %%

T_copv = T_copv_start;

rho_copv = double(py.CoolProp.CoolProp.PropsSI('D', 'P', P_copv, 'T', T_copv, 'Helium'));

m_he_copv = rho_copv * V_copv;
m_he_copv_start = m_he_copv;


T_ullage_lox = T_ullage_initial;
T_ullage_fuel = T_ullage_initial;

% we're assuming that the initial ullage is helium as well for simplicity
m_he_ullage_lox = P_ambient * V_ullage_lox / (R_He * T_ullage_lox);
m_he_ullage_fuel = P_ambient * V_ullage_fuel / (R_He * T_ullage_fuel);

U_ullage_lox = m_he_ullage_lox * cv_He * T_ullage_lox;
U_ullage_fuel = m_he_ullage_fuel * cv_He * T_ullage_fuel;

P_ullage_lox = P_ambient;
P_ullage_fuel = P_ambient;

%% INITIAL PRESSURIZATION (AMBIENT -> 535 psi)

m_he_ullage_lox_prefire = P_prefire * V_ullage_lox / (R_He * T_ullage_prefire);

m_he_ullage_fuel_prefire = P_prefire * V_ullage_fuel / (R_He * T_ullage_prefire);

% required helium for initial pressurization
m_he_prepress_lox = max(m_he_ullage_lox_prefire - m_he_ullage_lox, 0);
m_he_prepress_fuel = max(m_he_ullage_fuel_prefire - m_he_ullage_fuel, 0);
m_he_prepress = m_he_prepress_lox + m_he_prepress_fuel;

if m_he_copv - m_he_prepress <= 0
    error('COPV ran out of helium during initial pressurization.');
end

% Set the tank states directly to pre-fire condition.
m_he_ullage_lox = m_he_ullage_lox_prefire;
m_he_ullage_fuel = m_he_ullage_fuel_prefire;

T_ullage_lox = T_ullage_prefire;
T_ullage_fuel = T_ullage_prefire;

P_ullage_lox = P_prefire;
P_ullage_fuel = P_prefire;

U_ullage_lox = m_he_ullage_lox * cv_He * T_ullage_lox;
U_ullage_fuel = m_he_ullage_fuel * cv_He * T_ullage_fuel;

% Remove the prepress helium from the COPV and update its state using the
% same polytropic discharge model used during the burn.
rho_copv_old = rho_copv;
m_he_copv = m_he_copv - m_he_prepress;
rho_copv = m_he_copv / V_copv;

T_copv = T_copv * (rho_copv / rho_copv_old)^(polytropic_n - 1);

P_copv = double(py.CoolProp.CoolProp.PropsSI('P', 'D', rho_copv, 'T', T_copv, 'Helium'));

P_copv_ignition = P_copv;
T_copv_ignition = T_copv;
m_he_copv_ignition = m_he_copv;

fprintf('Assumed pre-fire ullage temp: %.1f K\n', T_ullage_prefire);
fprintf('Pre-fire LOX tank pressure:   %.1f psi\n', P_ullage_lox / PSI_TO_PA);
fprintf('Pre-fire fuel tank pressure:  %.1f psi\n', P_ullage_fuel / PSI_TO_PA);
fprintf('Helium used for prepress:      %.4f kg\n', m_he_prepress);
fprintf('COPV pressure at ignition:     %.1f psi\n', P_copv_ignition / PSI_TO_PA);
fprintf('COPV temperature at ignition:  %.1f K\n', T_copv_ignition);


%% BURN %%

N = ceil(burn_time / dt) + 1; % add buffer

time_array = zeros(1, N);
P_copv_array = zeros(1, N);
T_copv_array = zeros(1, N);
P_lox_array = zeros(1, N);
P_fuel_array = zeros(1, N);
T_ullage_lox_array = zeros(1, N);
T_ullage_fuel_array = zeros(1, N);
T_reg_out_array = zeros(1, N);
m_he_copv_array = zeros(1, N);
mdot_he_array = zeros(1, N);
rho_copv_array = zeros(1, N);
rho_ullage_lox_array = zeros(1,N);
rho_ullage_fuel_array = zeros(1,N);

% Initial burn values (immediately after initial press)
i = 1;
t = 0;

time_array(i) = t;
P_copv_array(i) = P_copv;
T_copv_array(i) = T_copv;
P_lox_array(i) = P_ullage_lox;
P_fuel_array(i) = P_ullage_fuel;
T_ullage_lox_array(i) = T_ullage_lox;
T_ullage_fuel_array(i) = T_ullage_fuel;
m_he_copv_array(i) = m_he_copv;
rho_copv_array(i) = rho_copv;

rho_ullage_lox_array(1) = m_he_ullage_lox / V_ullage_lox;
rho_ullage_fuel_array(1) = m_he_ullage_fuel / V_ullage_fuel;

m_he_burn_lox = 0;
m_he_burn_fuel = 0;
status = 0; % 0 = PASS, -1 = FAIL

while t < burn_time && status == 0
    dt_step = min(dt, burn_time - t);

    % Regulator outlet temperature at the loaded 480 psia condition.
    h_copv = double(py.CoolProp.CoolProp.PropsSI('H', 'P', P_copv, 'T', T_copv, 'Helium'));
    % h_in (copv) = h_out (reg outlet):
    T_reg_out = double(py.CoolProp.CoolProp.PropsSI('T', 'P', P_tank_burn, 'H', h_copv, 'Helium'));

    dm_lox_prop = min(mdot_lox * dt_step, m_lox);
    dm_fuel_prop = min(mdot_fuel * dt_step, m_fuel);

    V_ullage_lox_old = V_ullage_lox;
    V_ullage_fuel_old = V_ullage_fuel;

    m_lox = m_lox - dm_lox_prop;
    m_fuel = m_fuel - dm_fuel_prop;

    V_ullage_lox = V_lox_tank - m_lox / rho_lox;
    V_ullage_fuel = V_fuel_tank - m_fuel / rho_fuel;

    dV_lox = V_ullage_lox - V_ullage_lox_old;
    dV_fuel = V_ullage_fuel - V_ullage_fuel_old;

    dm_he_lox = 0;
    dm_he_fuel = 0;

    %% LOX ULLAGE
     if P_ullage_lox > P_tank_burn
        % The pre-fire 535 psi gas acts as a small buffer
        % No new helium is added yet
        T_ullage_lox = T_ullage_lox * (V_ullage_lox_old / V_ullage_lox)^(gamma_He - 1);
        P_ullage_lox = m_he_ullage_lox * R_He * T_ullage_lox / V_ullage_lox;

    else
        % Regulator is active. Hold the tank at 480 psi
        % Energy balance:
        % U_new = U_old + h_in*dm - P*dV
        % with:
        % U = m*cv*T
        % P*V = m*R*T
        % This gives an explicit value for the helium mass required

        U_old = m_he_ullage_lox * cv_He * T_ullage_lox;

        dm_he_lox = ((cv_He / R_He) * P_tank_burn * V_ullage_lox - U_old + P_tank_burn * dV_lox ) ...
            / (cp_He * T_reg_out);

        dm_he_lox = max(dm_he_lox, 0);
        m_he_ullage_lox = m_he_ullage_lox + dm_he_lox;

        T_ullage_lox = P_tank_burn * V_ullage_lox / (m_he_ullage_lox * R_He);

        P_ullage_lox = P_tank_burn;
    end

   
    %% FUEL ULLAGE %%
    if P_ullage_fuel > P_tank_burn
        % buffer active
        T_ullage_fuel = T_ullage_fuel * (V_ullage_fuel_old / V_ullage_fuel)^(gamma_He - 1);

        P_ullage_fuel = m_he_ullage_fuel * R_He * T_ullage_fuel / V_ullage_fuel;

    else
        % Regulator is active. Hold the tank at 480 psi
        U_old = m_he_ullage_fuel * cv_He * T_ullage_fuel;

        dm_he_fuel = ((cv_He / R_He) * P_tank_burn * V_ullage_fuel - U_old + P_tank_burn * dV_fuel ) ...
            / (cp_He * T_reg_out);

        dm_he_fuel = max(dm_he_fuel, 0);
        m_he_ullage_fuel = m_he_ullage_fuel + dm_he_fuel;

        T_ullage_fuel = P_tank_burn * V_ullage_fuel / ...
            (m_he_ullage_fuel * R_He);

        P_ullage_fuel = P_tank_burn;
    end

    % K factors
    

    % remove required helium from copv
    
    dm_he_total = dm_he_lox + dm_he_fuel;
    mdot_he = dm_he_total / dt_step;

    m_he_copv_new = m_he_copv - dm_he_total;

    if m_he_copv_new <= 0
        status = -1;
        fprintf('\nFAIL: COPV ran out of helium at t = %.3f s\n', t);
        break
    end

    rho_copv_old = rho_copv;
    m_he_copv = m_he_copv_new;
    rho_copv = m_he_copv / V_copv;

    T_copv = T_copv * (rho_copv / rho_copv_old)^(polytropic_n - 1);

    P_copv = double(py.CoolProp.CoolProp.PropsSI( ...
        'P', 'D', rho_copv, 'T', T_copv, 'Helium'));

    if P_copv <= P_tank_burn
        status = -1;
        fprintf('\nFAIL: COPV pressure is lower than tank pressure at t = %.3f s.\n', t + dt_step);
        break
    end

    m_he_burn_lox = m_he_burn_lox + dm_he_lox;
    m_he_burn_fuel = m_he_burn_fuel + dm_he_fuel;

    t = t + dt_step;
    i = i + 1;
    time_array(i) = t;
    P_copv_array(i) = P_copv;
    T_copv_array(i) = T_copv;
    P_lox_array(i) = P_ullage_lox;
    P_fuel_array(i) = P_ullage_fuel;
    T_ullage_lox_array(i) = T_ullage_lox;
    T_ullage_fuel_array(i) = T_ullage_fuel;
    T_reg_out_array(i) = T_reg_out;
    m_he_copv_array(i) = m_he_copv;
    mdot_he_array(i) = mdot_he;
    rho_copv_array(i) = rho_copv;
    rho_ullage_lox_array(i) = m_he_ullage_lox / V_ullage_lox;
    rho_ullage_fuel_array(i) = m_he_ullage_fuel / V_ullage_fuel;

end

time_array = time_array(1:i);
P_copv_array = P_copv_array(1:i);
T_copv_array = T_copv_array(1:i);
P_lox_array = P_lox_array(1:i);
P_fuel_array = P_fuel_array(1:i);
T_ullage_lox_array = T_ullage_lox_array(1:i);
T_ullage_fuel_array = T_ullage_fuel_array(1:i);
T_reg_out_array = T_reg_out_array(1:i);
m_he_copv_array = m_he_copv_array(1:i);
mdot_he_array = mdot_he_array(1:i);
rho_copv_array = rho_copv_array(1:i);
rho_ullage_lox_array = rho_ullage_lox_array(1:i);
rho_ullage_fuel_array = rho_ullage_fuel_array(1:i);

%% RESULTS & OUTPUT


m_he_burn_total = m_he_burn_lox + m_he_burn_fuel;
m_he_total_used = m_he_prepress + m_he_burn_total;

fprintf('\nCOPV Blowdown Results:\n');
fprintf('Burn time:                       %.2f s\n', burn_time);
fprintf('O/F:                             %.3f\n', OF);
fprintf('Total propellant mass:           %.2f kg\n', prop_mass);
fprintf('Total propellant mdot:           %.3f kg/s\n', mdot_total);
fprintf('Initial COPV Pressure:           %.1f psi\n', P_copv_start_psi);
fprintf('\n');
fprintf('Pre-fire tank pressure:          %.1f psi\n', P_prefire_psi);
fprintf('Steady-state tank pressure:      %.1f psi\n', P_tank_burn_psi);
fprintf('\n');
fprintf('COPV volume:                     %.2f L\n', V_copv_L);
fprintf('COPV pressure after ISO opens:   %.1f psi\n', P_copv_ignition / PSI_TO_PA);
fprintf('Final COPV pressure:             %.1f psi\n', P_copv / PSI_TO_PA);
fprintf('Final COPV temperature:          %.1f K\n', T_copv);
fprintf('\n');
fprintf('Initial COPV helium mass:        %.4f kg\n', m_he_copv_start);
fprintf('Helium used for prepress:        %.4f kg\n', m_he_prepress);
fprintf('Helium used during burn - LOX:   %.4f kg\n', m_he_burn_lox);
fprintf('Helium used during burn - fuel:  %.4f kg\n', m_he_burn_fuel);
fprintf('Total helium used:               %.4f kg\n', m_he_total_used);
fprintf('Helium remaining in COPV:        %.4f kg\n', m_he_copv);
fprintf('\n');



if status == 0
    fprintf('\nRESULT: PASS\n');
else
    fprintf('\nRESULT: FAIL\n');
end


%% PLOTS %%

figure;
plot(time_array, P_copv_array / PSI_TO_PA, 'LineWidth', 1.5);
hold on;
plot(time_array, P_lox_array / PSI_TO_PA, 'LineWidth', 1.5);
plot(time_array, P_fuel_array / PSI_TO_PA, '--', 'LineWidth', 1.5);
yline(P_tank_burn_psi, ':', '488 psia steady-state');
xlabel('Time [s]');
ylabel('Pressure [psi]');
title('Pressures During Burn');
legend('COPV', 'LOX Tank', 'Fuel Tank', 'Location', 'best');
grid on;

figure;
plot(time_array, T_copv_array, 'LineWidth', 1.5);
hold on;
plot(time_array, T_ullage_lox_array, 'LineWidth', 1.5);
plot(time_array, T_ullage_fuel_array, '--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Temperature [K]');
title('Gas Temperatures During Burn');
legend('COPV', 'LOX Ullage', 'Fuel Ullage', 'Location', 'best');
grid on;



figure;

yyaxis left
plot(time_array, rho_copv_array, 'LineWidth', 1.5);
ylabel('COPV Density [kg/m^3]');

yyaxis right
plot(time_array, rho_ullage_lox_array, 'LineWidth', 1.5);
hold on;
plot(time_array, rho_ullage_fuel_array, '--', 'LineWidth', 1.5);
ylabel('Ullage Density [kg/m^3]');

xlabel('Time [s]');
title('Helium Density During Burn');
legend('COPV', 'LOX Ullage', 'Fuel Ullage', 'Location', 'best');
grid on;