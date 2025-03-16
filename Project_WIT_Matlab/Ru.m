function value = Ru(Cu, Cv, V_mu, K_mu, K_mv)
    value = (V_mu*Cu)./((K_mu+Cu).*(1+Cv/K_mv));
    %disp(value)
end
