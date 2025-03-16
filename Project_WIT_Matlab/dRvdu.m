function value = dRvdu(cu, cv, V_mu, V_mfv, K_mu, K_mv, K_mfu,rq)
    value = rq*dRudu(cu, cv, V_mu, K_mu, K_mv)-diag(V_mfv*(1./(K_mfu*(1+cu/K_mfu).^2)));
end
