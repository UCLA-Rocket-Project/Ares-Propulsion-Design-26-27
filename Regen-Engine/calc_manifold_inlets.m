load('Geo_struct.mat', 'Geo');

%% Aft Manifold Geometry
theta = linspace(0, pi, 181); % Row vector
R_boattail = 0.0635; % m, inlet right above exit
%R_boattail = 0.06530975; % m, inlet 2mm above exit
R_outer_exit = Geo.r_outer_wall(end-1);
Z_exit = Geo.pos_i(end-1);

D_downcomer = convlength(0.625, 'in', 'm'); % m, -12 AN
A_downcomer = pi * (D_downcomer/2)^2;
A_manifold_inlet = A_downcomer; % max
A_manifold_min = 2*(Geo.w_channel(end-1) * Geo.h_channel(end-1)); 

W_max = 0.85*(R_boattail - R_outer_exit); % knockdown factor to make downcomer fit more within boattail and 1mm raised height
W_min = Geo.w_channel(end-1); 

% Calculate required Area and Width
A_theta = A_manifold_inlet .* (1 - (theta ./ pi)) + A_manifold_min .* (theta ./ pi);
W_theta = W_max .* (1 - (theta ./ pi)) + W_min .* (theta ./ pi);
writematrix(A_theta.', 'Manifold_Area.csv');
writematrix(W_theta.', 'Manifold_Width.csv');

%% Inlet Geometry (in mm)
A_downcomer_mm = A_downcomer * 1e6;
A_inlet = A_downcomer_mm / (0.5 * Geo.num_channel);
W_inlet = 3.83998387; % from CAD
proj_angle = convang(73, 'deg', 'rad'); % relative nozzle exit slant rel. to horizontal
roof_angle = convang(45, 'deg', 'rad'); % angle for pointed top (prevent printing overhang)
r_fillet = 0.5; % CAD fillet radius
r_fillet_top = 1; % to increase rectangular feature more
A_loss_rs = r_fillet^2 * (1 + pi/4); % Area loss to right angle fillets (small radius)
A_loss_rl = r_fillet_top^2 * (1 + pi/4);
A_loss_o = (r_fillet/2)^2 * (cot((roof_angle+pi/2)/2) - (pi - (roof_angle + pi/2))/2); % Area loss to obtuse angle fillets
A_projected = A_inlet * sin(proj_angle);
A_theo = A_projected + 2*A_loss_rs + A_loss_rl + 2*A_loss_o;
h_tri = 0.5 * W_inlet * tan(roof_angle);
A_tri = 0.5 * W_inlet * h_tri;
A_rec = A_theo - A_tri;
h_rec = A_rec / W_inlet; % input in CAD



