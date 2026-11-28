function Stress = stressAnalysis(Geo, Loop, Mat, Param, d, Temp)
    p_diff = Loop.P_loc - Param.Pc;
    r_in = Loop.D_g_loc/2;
    r_out = r_in+Geo.wall_thickness(d)+Loop.ch+Geo.out_wall_thickness;
    r_avg = (r_in + r_out)/2;
    Mat.E = interp1(Mat.ref_temps, Mat.E_ref, Temp.T_cw, 'linear', 'extrap');
    Mat.nu = interp1(Mat.ref_temps, Mat.nu_ref, Temp.T_cw, 'linear', 'extrap');
    Mat.alpha = interp1(Mat.alpha_ref_temps, Mat.alpha_ref, Temp.T_cw, 'linear', 'extrap');
    
    % Shear
    Stress.sigma_s_mech = 0.5*p_diff*(Loop.cw/Geo.wall_thickness(d));
    Stress.sigma_s_therm = ((Loop.cw+Geo.w_rib)*Mat.alpha*Mat.E*...
        (Loop.T_iw-Loop.T_ow)*Geo.wall_thickness(d)*Geo.out_wall_thickness)/...
        (5*Geo.w_rib*(1-Mat.nu)*...
        (Geo.wall_thickness(d)+Geo.out_wall_thickness)^2);
    Stress.sigma_s_tot = Stress.sigma_s_mech + Stress.sigma_s_therm;
    
    % Hoop
    Stress.sigma_h_mech = (p_diff * r_avg)/...
        (Geo.wall_thickness(d)+Geo.out_wall_thickness);
    Stress.sigma_h_therm = (Mat.alpha*Mat.E*(Loop.T_iw-Loop.T_ow)*...
        Geo.out_wall_thickness) / ((1-Mat.nu)*...
        (Geo.wall_thickness(d) + Geo.out_wall_thickness));
    Stress.sigma_h_tot = Stress.sigma_h_mech + Stress.sigma_h_therm;
    
    % Bending
    Stress.sigma_b = (p_diff * Loop.cw^2)/(2*Geo.wall_thickness(d)^2);
    
    % Axial
    if (d < Geo.pos_throat)
        Stress.sigma_a = (Param.Pc*pi*((Loop.D_g_loc/2)^2-Geo.Rt^2))/...
            ((pi*(r_out^2 - Geo.Rt^2))-...
            (Geo.num_channel*Loop.cw*Loop.ch));
    else
        Stress.sigma_a = 0;
    end
    
    Stress.sigma_VM = sqrt(0.5*((Stress.sigma_a-Stress.sigma_h_tot)^2+...
        (Stress.sigma_h_tot - Stress.sigma_b)^2 +...
        (Stress.sigma_b - Stress.sigma_a)^2) + 3*Stress.sigma_s_tot^2);
end