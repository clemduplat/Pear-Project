function value = dRudu(cu, cv, V_mu, K_mu, K_mv)
    vals = (V_mu*K_mu)./((1. + cv/K_mv).*(K_mu + cu).^2);
    value = diag(vals);
end