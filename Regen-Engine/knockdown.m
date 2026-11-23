clc; clear; close all;

load('Geo_struct.mat');
T_measured = 394.9; % from TC data, bulk temp (K)
Loc_measured = 0.32; % m, from the exit
Loc_script = Geo.pos_i(end) - Loc_measured; % m, relative to nozzle exit
[~, node_script] = min(abs(Geo.pos_i - Loc_script)); % find closest axial station to measured loc


tol_kd = 1e-1; % K
knockdown_iter = 0;
kd_prev_error = 0;
kd_guess = 1;
kd_prev_guess = 0.8;
while knockdown_iter < 50
    knockdown_iter = knockdown_iter + 1;
    Arrays = regen_back_calc(kd_guess, false);
    T_calc = Arrays.T_bulk_array(node_script);

    kd_error = T_calc - T_measured;
    fprintf('Iter %d: kd_guess = %.4f | T_bulk = %.2f K | T_hw = %.2f K | T_cw = %.2f K | Error = %.3f K\n', ...
        knockdown_iter, kd_guess, T_calc, Arrays.T_hw_array(node_script), Arrays.T_cw_array(node_script), kd_error);
    if abs(kd_error) < tol_kd
        fprintf('Converged on knockdown factor: %.4f\n', kd_guess);
        break;
    end
    if knockdown_iter == 1
        kd_prev_error = kd_error;
        temp_kd = kd_guess;
        kd_guess = kd_prev_guess;
        kd_prev_guess = temp_kd;
    else
        kd_next = kd_guess - kd_error * (kd_guess - kd_prev_guess) / (kd_error - kd_prev_error + 1e-10);
        kd_next = max(0.05, min(kd_next, 1.2)); % bound the secant guess to reasonable values
        kd_prev_guess = kd_guess;
        kd_prev_error = kd_error;
        kd_guess = kd_next;
    end
end

Arrays = regen_back_calc(kd_guess, true);
