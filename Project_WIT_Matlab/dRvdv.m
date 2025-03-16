function value = dRvdv(cu, cv, V_mu, K_mu, K_mv,rq)
    value = rq*dRudv(cu,cv,V_mu,K_mu,K_mv);
end