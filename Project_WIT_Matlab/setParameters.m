function return_list = setParameters(storage, p_atm, R_g, Ea_vmfv_ref, Ea_vmu_ref, T_ref, V_mu_ref, V_mfv_ref)
    switch storage
        case 1
            T_cel = 25; eta_u = 20.8; eta_v = 0.04;
        case 2
            T_cel = 20; eta_u = 20.8; eta_v = 0;
        case 3
            T_cel = 7; eta_u = 20.8; eta_v = 0;
        case 4
            T_cel = -1; eta_u = 20.8; eta_v = 0;
        case 5
            T_cel = -1; eta_u = 2; eta_v = 5;
        case 6
            T_cel = -1; eta_u = 2; eta_v = 0.7;
        otherwise
            error('Error! No valid storage conditions are specified');
    end

    T = T_cel + 273.15;
    Cu_amb = (p_atm * eta_u / 100) / (R_g * T);
    Cv_amb = (p_atm * eta_v / 100) / (R_g * T);

    V_mu = V_mu_ref * exp(Ea_vmu_ref / R_g * (1 / T_ref - 1 / T));
    V_mfv = V_mfv_ref * exp(Ea_vmfv_ref / R_g * (1 / T_ref - 1 / T));

    return_list = [T_cel, eta_u, eta_v, T, Cu_amb, Cv_amb, V_mu, V_mfv];
    
end
