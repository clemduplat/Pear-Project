function value = Rv(Cu, Cv, V_mu, K_mu, K_mv, K_mfu, rq,V_mfv)
    value = rq*Ru(Cu, Cv, V_mu, K_mu, K_mv) + V_mfv./(1 + Cu./K_mfu);
end