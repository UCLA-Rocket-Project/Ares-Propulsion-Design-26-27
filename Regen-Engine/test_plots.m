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