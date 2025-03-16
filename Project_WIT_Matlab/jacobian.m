function matrix = jacobian(Ku, Kv, Fu, cu, cv, V_mu, K_mu, K_mv, Hu1, Fv1, K_mfu, Hv1, V_mfv,rq)
    r1Du = Ku - Fu*dRudu(cu, cv, V_mu, K_mu, K_mv) + Hu1;
    r1Dv = - Fu*dRudv(cu, cv, V_mu, K_mu, K_mv);
    r2Du = - Fv1*dRvdu(cu, cv, V_mu, V_mfv, K_mu, K_mv, K_mfu,rq);
    r2Dv = Kv - Fv1*dRvdv(cu, cv, V_mu, K_mu, K_mv,rq) + Hv1;
    matrix = [[(r1Du) (r1Dv)]; [(r2Du) (r2Dv)]];
    %size(matrix)
end
