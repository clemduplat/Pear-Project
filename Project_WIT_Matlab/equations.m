function [rows] = equations(cu, cv, Ku, Kv, Fu, Fv1, Hu1, Hu2, Hv1, Hv2, V_mu, K_mu, K_mv, K_mfu, rq, V_mfv)
    row1 = Ku*cu - Fu*Ru(cu, cv, V_mu, K_mu, K_mv) + Hu1*cu - Hu2;
    row2 = Kv*cv - Fv1*Rv(cu, cv, V_mu, K_mu, K_mv, K_mfu, rq, V_mfv) + Hv1*cv - Hv2;
    %disp(Ru(cu, cv, V_mu, K_mu, K_mv))
    %disp(Rv(cu, cv, V_mu, K_mu, K_mv, K_mfu, rq, V_mfv))
    rows = [row1; row2];
    %size(rows)
end
