function Cool = load_coolprop(table_filename, table_version_req, x_eth, x_h2o, mfrac_eth, mfrac_h2o)
    % Define Table bounds (Tmax < Loop.T_sat at Pmin or coolprop will crash)
    tables_loaded = false;
    if isfile(table_filename)
        tbl = load(table_filename);
        if isfield(tbl, 'table_version') && tbl.table_version == table_version_req
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
end