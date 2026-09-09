function Temp = temp_iteration(Param, Cantera, Y_str, Geo, Gas, Cool, Mat, Loop, d) % HW Temp Iteration
    Temp.T_tc_guess = 1500; % Initial Guess (1.25 FOS applied to material melting point)
    T_tc_prev_guess = 1000; % Cooler lower bound (Luca's value)
    tol_q = 0.1; % w
    q_prev_error = 0;
    Temp.iter_T = 0;
    while true
        Temp.iter_T = Temp.iter_T + 1;

        % Gas convection HTC with Bartz
        Temp.sigma = 1 / ...
            ((0.5 * (Temp.T_tc_guess/Gas.T_stag)*(1 + (Gas.gamma-1)/2 * Gas.M_local(d)^2) + 0.5)^0.68 *...
            (1 + (Gas.gamma-1)/2 * Gas.M_local(d)^2)^0.12);

        Temp.h_i = ((0.026/Geo.D_t^0.2)*...
            (Gas.mu_g_local(d)^0.2/Gas.prandtl_g_local(d)^0.6)*...
            (Param.Pc/Param.cstar_act)^0.8)*...
            (Geo.D_t/Geo.R_curve)^0.1*...
            (Geo.At/Gas.A_local(d))^0.9*...
            Temp.sigma;
        
        i_w = get_iw(Cantera, Y_str, Gas, Temp);
        Temp.q_eq = Temp.h_i * Loop.A_g_loc * (Gas.i_aw(d) - i_w);

        Temp.k_Al2O3_loc = interp1(Mat.k_Al2O3_ref_temps, Mat.k_Al2O3_ref, Temp.T_tc_guess, 'linear', 'extrap');
        Temp.k_ZrO2_loc = interp1(Mat.k_ZrO2_ref_temps, Mat.k_ZrO2_ref, Temp.T_tc_guess, 'linear', 'extrap');
        Temp.k_coating = 0.9722*Temp.k_ZrO2_loc + 0.02778*Temp.k_Al2O3_loc; % Volume Fraction
        
        R_tc = (log((Loop.D_g_loc/2 + Geo.coat_thickness)/(Loop.D_g_loc/2)))/...
            (2*pi*Temp.k_coating*Geo.dl(d));
        Temp.T_hw = Temp.T_tc_guess - Temp.q_eq*R_tc;
        Temp.T_hw = max(Loop.T_bulk + 1, min(Temp.T_hw, Loop.Taw_loc - 1));

        % Fin Efficiency
        Temp.k_w_loc = interp1(Mat.k_w_ref_temps, Mat.k_w_ref, Temp.T_hw, 'linear', 'extrap');
        Temp.fin_m = sqrt((2*Loop.h_c)/(Temp.k_w_loc * Geo.w_rib));
        Temp.fin_eff = tanh(Temp.fin_m * Loop.ch) / (Temp.fin_m * Loop.ch);
        Temp.h_c_f = Loop.h_c*(Loop.cw+2*Temp.fin_eff*Loop.ch)/(Loop.cw+Geo.w_rib); % Fin corrected Loop.h_c
        
        Temp.T_cw = Temp.T_hw - (Temp.q_eq)*...
            log(1+2*Geo.wall_thickness(d)/(Loop.D_g_loc+2*Geo.coat_thickness))/(2*pi*Geo.dl(d)*Temp.k_w_loc); % Cold wall temp derived from guess

        % coolant side shares the convection term in both regimes so q_c is continuous at t_cw = t_sat
        mass_flux = Param.mdot_f / Loop.A_conv;
        dT_sub = Loop.T_sat - Loop.T_bulk;
        dT_bulk = Temp.T_cw - Loop.T_bulk;
        dT_sat = Temp.T_cw - Loop.T_sat;
        Stanton = abs(Temp.q_eq/Loop.A_base_loc) / (mass_flux * Loop.cp_c * dT_sub); % dimensionless
        Peclet = mass_flux * Loop.D_h_loc * Loop.cp_c / Loop.k_c;
        Boiling = abs(Temp.q_eq/Loop.A_base_loc) / (mass_flux * Loop.h_fg);

        Temp.h_nb = 0;
        Temp.h_tp = Temp.h_c_f;
        Temp.h_c_FEA = Loop.h_c;
        % Zhu-Bi-Yan Correlation
        q_c = Temp.h_c_f*Loop.A_base_loc*(Temp.T_cw - Loop.T_bulk);
        if (Temp.T_cw >= Loop.T_sat)
            if Stanton <= 38*Peclet^(-0.38) % partially boiling
                Temp.h_nb = 0.00122*...
                    (((Loop.k_c^0.79)*(Loop.cp_c^0.45)*(Loop.rho_c_l^0.49))/...
                    ((Loop.surften^0.5)*(Loop.mu_c^0.29)*(Loop.h_fg^0.24)*(Loop.rho_c_v^0.24)))*...
                    ((dT_sat)^0.24)*...
                    (max(0, Cool.get_P_sat(min(Temp.T_cw, 513)) - Loop.P_loc))^0.75; % Cap at critical pressure for ethanol
                S = (1 / (1 + 0.055*Loop.Re^0.16)) * (Loop.T_sat/dT_sub)^0.28;
                Temp.h_tp = sqrt(Temp.h_c_f^2 + (S*Temp.h_nb*dT_sat/dT_bulk)^2);
                Temp.h_c_FEA = sqrt(Loop.h_c^2 + (S*Temp.h_nb*dT_sat/dT_bulk)^2); % without fin efficiency, for FEA
                q_c = Temp.h_tp * Loop.A_base_loc*(Temp.T_cw - Loop.T_bulk);
            else % fully developed boiling
                q_c = 1000 * Loop.A_base_loc * dT_sat / (32 * exp(-Loop.P_loc/8.6*10^6)); % W, Zhu-Bi
            end
        end

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
            temp_T = Temp.T_tc_guess;
            Temp.T_tc_guess = T_tc_prev_guess;
            T_tc_prev_guess = temp_T;
        else
            T_next = Temp.T_tc_guess - Temp.q_error * (Temp.T_tc_guess - T_tc_prev_guess) / (Temp.q_error - q_prev_error + 1e-10); % small buffer to avoid divide by 0
            T_next = max(Loop.T_bulk + 1, min(T_next, Loop.Taw_loc - 1));
            T_tc_prev_guess = Temp.T_tc_guess;
            q_prev_error = Temp.q_error;
            Temp.T_tc_guess = T_next;
        end
        
    end
    Temp.h_g = Temp.q_eq / (Loop.A_g_loc * (Gas.Taw(d) - Temp.T_tc_guess)); % for FEA

    function i_w = get_iw(Cantera,Y_str, Gas, Temp)
        Cantera.TPY = py.tuple({Temp.T_tc_guess, Gas.pressure(d), Y_str});
        Cantera.equilibrate('TP');
        i_w = double(Cantera.enthalpy_mass);
    end
end