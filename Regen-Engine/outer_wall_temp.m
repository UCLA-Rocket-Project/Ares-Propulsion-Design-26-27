function Temp = outer_wall_temp(Geo, Loop, Temp, Param, d)
    R_ow = Loop.D_g_loc/2 + Geo.wall_thickness(d) + Loop.ch + Geo.out_wall_thickness;
    C = (Temp.k_w_loc*Geo.w_rib*Temp.fin_m)/sinh(Temp.fin_m*Loop.ch);
    Temp.T_fin_tip = (C*Temp.T_cw - C*Loop.T_bulk + C*Loop.T_bulk*cosh(Temp.fin_m*Loop.ch) +...
        Loop.h_c*Loop.cw*Loop.T_bulk)/(C*cosh(Temp.fin_m*Loop.ch)+Loop.h_c*Loop.cw);
    
    T_ow_error = realmax;
    tol = 0.1;
    T_ow_guess = Temp.T_fin_tip;
    while (abs(T_ow_error) > tol)
        k_visc = 15.51 * 10^-6; % for air at 298K and 1 atm
        alpha = 22.39 * 10^-6;
        Pr = k_visc/alpha;
        beta = 3.43*10^-3;
        Gr = (9.81*beta*(T_ow_guess - Param.T_amb)*(2*R_ow)^3)/k_visc^2;
        Ra = Gr * Pr;
        Nu = (0.6 + (0.387*Ra^(1/6))/((1+(0.559/Pr)^(9/16))^(8/27)))^2; % Churchill & Chu
        k = 26.23*10^-3;
        h_a = (Nu*k)/(R_ow*2);
        R_cond = Geo.out_wall_thickness/(Temp.k_w_loc*(Geo.w_rib+Loop.cw));
        R_conv = 1/(h_a*(Geo.w_rib+Loop.cw));
        Temp.T_ow = Param.T_amb + (R_conv*(Temp.T_fin_tip-Param.T_amb))/(R_cond+R_conv);
    
        T_ow_error = Temp.T_ow - T_ow_guess;
        T_ow_guess = T_ow_guess + 0.5*T_ow_error;
    end
    Temp.T_ow = T_ow_guess;
end