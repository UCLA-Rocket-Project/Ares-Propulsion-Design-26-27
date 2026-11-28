% Elfaham Correlation
% Find Ms correction term
% if 1e-10 <= Boiling && Boiling < 1e-3 % actual lower bound 1e-5, change to this if Boiling within appropriate range
%     Ms = 0.7;
% elseif 1e-3 <= Boiling && Boiling < 5e-3
%     Ms = 1.5;
% elseif 5e-3 <= Boiling && Boiling < 1e-2
%     Ms = 1.3;
% else
%     Ms = 1.1;
% end
% q_c = Temp.h_c_f*Loop.A_base_loc*(Temp.T_cw - Loop.T_bulk);
% Temp.h_nb = 0;
% if Temp.T_cw > Loop.T_sat
%     Temp.h_nb = 55 * (abs(Temp.q_eq/Loop.A_base_loc))^0.67 * Loop.P_reduced^0.12 * (-log10(Loop.P_reduced))^(-0.55) * Param.MW_coolant^(-0.5);
%     E = (1 + Loop.quality*Loop.prandtl_c*(Loop.rho_c_l/Loop.rho_c_v - 1))^0.35;
%     S = Ms / (1 + 0.055 * E^0.1 * Loop.Re^0.16);
%     Temp.h_tp = sqrt((S*Temp.h_nb)^2 + (E*Temp.h_c_f)^2);
%     q_c = Temp.h_tp * Loop.A_base_loc*(Temp.T_cw - Loop.T_bulk);
% end

% Chen Correlation
% q_c = Temp.h_c_f*Loop.A_base_loc*(Temp.T_cw - Loop.T_bulk);
% if Temp.T_cw > Loop.T_sat
%     Temp.h_nb = 0.00122*...
%         (((Loop.k_c^0.79)*(Loop.cp_c^0.45)*(Loop.rho_c_l^0.49))/...
%         ((Loop.surften^0.5)*(Loop.mu_c^0.29)*(Loop.h_fg^0.24)*(Loop.rho_c_v^0.24)))*...
%         ((Temp.T_cw-Loop.T_sat)^0.24)*...
%         (max(0, Cool.get_P_sat(min(Temp.T_cw, 513)) - Loop.P_loc))^0.75; % Cap at critical pressure for ethanol
%     S = 1/(1+(2.53*(10^-6))*(Loop.Re^1.17));
%     q_c = q_c + S*Temp.h_nb*(Temp.h_c_f/Loop.h_c)*Loop.A_base_loc*(Temp.T_cw - Loop.T_sat);
% end