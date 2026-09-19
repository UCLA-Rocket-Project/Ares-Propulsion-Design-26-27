%% Aft Manifold Geometry
r_channel_base = Geo.D_channel_base ./ 2;
r_outer_jacket = r_channel_base + Geo.h_channel;
r_outer_wall = r_channel_base + Geo.h_channel + Geo.out_wall_thickness;

theta = linspace(0, pi, 100); % Row vector
D_boattail = 0.13025; % m
R_boattail = D_boattail / 2;
R_outer_exit = r_outer_wall(end);
Z_exit = Geo.pos_i(end);

% FIX: Ensure this is the radius, not the diameter, or the area inflates by 400%
A_downcomer = pi * (0.01332)^2; % Assuming 26.64 mm was diameter -> 13.32 mm radius
A_manifold_inlet = A_downcomer / 2; 
A_manifold_min = Geo.w_channel(end) * Geo.h_channel(end); 

W_max = R_boattail - R_outer_exit;
W_min = Geo.w_channel(end); 

% 1. Linear schedules for required Area and Width
A_theta = A_manifold_inlet .* (1 - (theta ./ pi)) + A_manifold_min .* (theta ./ pi);
W_theta = W_max .* (1 - (theta ./ pi)) + W_min .* (theta ./ pi);

% 2. Iterative Integral Solver for Conformal Height
Z_top_array = zeros(size(theta));
R_top_array = zeros(size(theta));

for k = 1:length(theta)
    A_target = A_theta(k);
    W_loc = W_theta(k);
    R_out_base = R_outer_exit + W_loc; % The boattail rim at this theta
    
    A_calc = 0;
    % March backwards from the exit plane up the nozzle
    for idx = length(Geo.pos_i)-1 : -1 : 1
        z_test = Geo.pos_i(idx:end);
        r_in = r_outer_wall(idx:end);
        
        % 1. Normalized axial distance (0 at top of manifold, 1 at bottom exit plane)
        norm_z = (z_test - z_test(1)) ./ (z_test(end) - z_test(1));
        
        % 2. Straight line base (what you previously had)
        r_straight = r_in(1) + (R_out_base - r_in(1)) .* norm_z;
        
        % 3. Add a parabolic bulge (convex outward)
        bulge_factor = 0.6; % 0 = flat roof. 1.0 = highly curved teardrop.
        r_out_curve = r_straight + bulge_factor .* (R_out_base - r_in(1)) .* norm_z .* (1 - norm_z);
        
        % 4. Calculate exact cross-sectional area between curved roof and nozzle
        A_calc = trapz(z_test, r_out_curve - r_in);
        
        if A_calc >= A_target
            Z_top_array(k) = z_test(1);
            R_top_array(k) = r_in(1);
            break;
        end
    end
    
    % Fallback if geometry forces an impossible area
    if A_calc < A_target
        Z_top_array(k) = Geo.pos_i(1);
        R_top_array(k) = r_outer_wall(1);
    end
end

% 3. 3D Cartesian Coordinates for CAD Guide Curves
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

% Conformal Plotting - Inlet (0 degrees)
idx_0 = find(Geo.pos_i >= Z_top_array(1), 1);
z_curve_0 = Geo.pos_i(idx_0:end);
r_curve_0 = r_outer_wall(idx_0:end);
% Draw sequence: Top point -> Sloped Roof -> Bottom Floor -> Curved Inner Wall
Z_man_0 = [Z_top_array(1), Z_exit, z_curve_0(end:-1:1)];
R_man_0 = [R_top_array(1), R_outer_exit + W_theta(1), r_curve_0(end:-1:1)];
plot(Z_man_0, R_man_0, 'm-', 'LineWidth', 2, 'DisplayName', 'Manifold Inlet (0^\circ)');

% Conformal Plotting - End (180 degrees)
idx_180 = find(Geo.pos_i >= Z_top_array(end), 1);
z_curve_180 = Geo.pos_i(idx_180:end);
r_curve_180 = -r_outer_wall(idx_180:end); % Negative for bottom half of plot
Z_man_180 = [Z_top_array(end), Z_exit, z_curve_180(end:-1:1)];
R_man_180 = [-R_top_array(end), -(R_outer_exit + W_theta(end)), r_curve_180(end:-1:1)];
plot(Z_man_180, R_man_180, 'm--', 'LineWidth', 2, 'DisplayName', 'Manifold End (180^\circ)');

% Manifold Top Edge Projection (Connecting 0 to 180 degrees)
plot(Z_top_array, R_top_array .* cos(theta), 'm-', 'LineWidth', 2, 'DisplayName', 'Manifold Top Edge Profile');

% Exit Plane Vertical Line (Outer wall cap)
plot([Z_exit, Z_exit], [R_outer_exit, -R_outer_exit], 'm--', 'LineWidth', 2, 'HandleVisibility', 'off');

title('Regen 1D Profile')
xlabel('Axial Position x (m)');
ylabel('Radial Position y (m)');
axis equal;
xline(0, 'r--', 'Throat', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
exportgraphics(gcf, 'geometry.pdf', 'ContentType','vector');
hold off;

%% Aft Manifold Geometry
r_channel_base = Geo.D_channel_base ./ 2;
r_outer_jacket = r_channel_base + Geo.h_channel;
r_outer_wall = r_channel_base + Geo.h_channel + Geo.out_wall_thickness;

theta = linspace(0, pi, 100); % Row vector
D_boattail = 0.13025; % m
R_boattail = D_boattail / 2;
R_outer_exit = r_outer_wall(end);
Z_exit = Geo.pos_i(end);

% Fixed to use Radius, preventing 400% area inflation
A_downcomer = pi * (0.01332)^2; 
A_manifold_inlet = A_downcomer / 2; 
A_manifold_min = Geo.w_channel(end) * Geo.h_channel(end); 

W_max = R_boattail - R_outer_exit;
W_min = Geo.w_channel(end); 

% 1. Linear schedules for required Area and Width
A_theta = A_manifold_inlet .* (1 - (theta ./ pi)) + A_manifold_min .* (theta ./ pi);
W_theta = W_max .* (1 - (theta ./ pi)) + W_min .* (theta ./ pi);

% 2. Iterative Integral Solver for Conformal Height
Z_top_array = zeros(size(theta));
R_top_array = zeros(size(theta));

bulge_factor = 0.6; % 0 = flat sloped roof, >0 = convex teardrop roof

for k = 1:length(theta)
    A_target = A_theta(k);
    W_loc = W_theta(k);
    R_out_base = R_outer_exit + W_loc; 
    
    A_prev = 0;
    % March backwards from the exit plane up the nozzle
    for idx = length(Geo.pos_i)-1 : -1 : 1
        z_test = Geo.pos_i(idx:end);
        r_in = r_outer_wall(idx:end);
        
        % Protect against division by zero at the exact exit plane
        if (z_test(end) - z_test(1)) == 0
            norm_z = zeros(size(z_test));
        else
            norm_z = (z_test - z_test(1)) ./ (z_test(end) - z_test(1));
        end
        
        r_straight = r_in(1) + (R_out_base - r_in(1)) .* norm_z;
        r_out_curve = r_straight + bulge_factor .* (R_out_base - r_in(1)) .* norm_z .* (1 - norm_z);
        
        A_calc = trapz(z_test, r_out_curve - r_in);
        
        if A_calc >= A_target
            % Sub-grid interpolation for a mathematically smooth curve
            if A_calc == A_prev
                frac = 0;
            else
                frac = (A_target - A_prev) / (A_calc - A_prev);
            end
            
            Z_prev = Geo.pos_i(idx+1);
            R_prev = r_outer_wall(idx+1);
            Z_curr = Geo.pos_i(idx);
            R_curr = r_outer_wall(idx);
            
            Z_top_array(k) = Z_prev + frac * (Z_curr - Z_prev);
            R_top_array(k) = R_prev + frac * (R_curr - R_prev);
            break;
        end
        A_prev = A_calc;
    end
    
    if A_calc < A_target
        Z_top_array(k) = Geo.pos_i(1);
        R_top_array(k) = r_outer_wall(1);
    end
end

blend_percent = 0.80; 
blend_idx = floor(blend_percent * length(theta));

Z_blend_start = Z_top_array(blend_idx);
R_blend_start = R_top_array(blend_idx);

% Save the true solver endpoints for the 180-degree channel!
Z_end_target = Z_top_array(end);
R_end_target = R_top_array(end);

% Normalized blend parameter (0 at blend start, 1 at dead end)
t_blend = linspace(0, 1, length(theta) - blend_idx + 1);

% Quadratic smoothing targeting the true final channel geometry, NOT zero
Z_smooth = Z_end_target - (Z_end_target - Z_blend_start) .* (1 - t_blend).^2;
R_smooth = R_end_target - (R_end_target - R_blend_start) .* (1 - t_blend).^2;

Z_top_array(blend_idx:end) = Z_smooth;
R_top_array(blend_idx:end) = R_smooth;

% 3. 3D Cartesian Coordinates for CAD Guide Curves
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

% Conformal Plotting - Inlet (0 degrees)
W_0 = W_theta(1);
R_out_base_0 = R_outer_exit + W_0;
z_roof_0 = linspace(Z_top_array(1), Z_exit, 20);
if (z_roof_0(end) - z_roof_0(1)) == 0; norm_z_0 = zeros(size(z_roof_0)); else; norm_z_0 = (z_roof_0 - z_roof_0(1)) ./ (z_roof_0(end) - z_roof_0(1)); end
r_straight_0 = R_top_array(1) + (R_out_base_0 - R_top_array(1)) .* norm_z_0;
r_roof_0 = r_straight_0 + bulge_factor .* (R_out_base_0 - R_top_array(1)) .* norm_z_0 .* (1 - norm_z_0);

idx_0 = find(Geo.pos_i >= Z_top_array(1), 1);
z_inner_0 = Geo.pos_i(idx_0:end);
r_inner_0 = r_outer_wall(idx_0:end);
if z_inner_0(1) > Z_top_array(1)
    z_inner_0 = [Z_top_array(1), z_inner_0];
    r_inner_0 = [R_top_array(1), r_inner_0];
end
Z_man_0 = [z_roof_0, z_inner_0(end:-1:1)];
R_man_0 = [r_roof_0, r_inner_0(end:-1:1)];
plot(Z_man_0, R_man_0, 'm-', 'LineWidth', 2, 'DisplayName', 'Manifold Inlet (0^\circ)');

% Conformal Plotting - End (180 degrees)
W_180 = W_theta(end);
R_out_base_180 = R_outer_exit + W_180;
z_roof_180 = linspace(Z_top_array(end), Z_exit, 20);
if (z_roof_180(end) - z_roof_180(1)) == 0; norm_z_180 = zeros(size(z_roof_180)); else; norm_z_180 = (z_roof_180 - z_roof_180(1)) ./ (z_roof_180(end) - z_roof_180(1)); end
r_straight_180 = R_top_array(end) + (R_out_base_180 - R_top_array(end)) .* norm_z_180;
r_roof_180 = r_straight_180 + bulge_factor .* (R_out_base_180 - R_top_array(end)) .* norm_z_180 .* (1 - norm_z_180);

idx_180 = find(Geo.pos_i >= Z_top_array(end), 1);
z_inner_180 = Geo.pos_i(idx_180:end);
r_inner_180 = r_outer_wall(idx_180:end);
if z_inner_180(1) > Z_top_array(end)
    z_inner_180 = [Z_top_array(end), z_inner_180];
    r_inner_180 = [R_top_array(end), r_inner_180];
end
Z_man_180 = [z_roof_180, z_inner_180(end:-1:1)];
R_man_180 = [-r_roof_180, -r_inner_180(end:-1:1)]; % Negative for bottom half
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