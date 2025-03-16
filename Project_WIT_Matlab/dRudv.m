function value = dRudv(cu, cv, V_mu, K_mu, K_mv)
    vals = -(K_mv*V_mu*cu)./((K_mu+cu).*(K_mv+cv).^2);
    value = diag(vals);
end
