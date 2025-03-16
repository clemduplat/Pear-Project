function [Cu, Cv] = newton_raphson(Cu0, Cv0, tolerance, max_iterations, Ku, Kv, Fu, V_mu, K_mu, K_mv, Hu1, Hu2, Fv1, K_mfu, Hv1, Hv2, V_mfv, rq)
    res = tolerance+1;
    Cu = Cu0;
    Cv = Cv0;
    C = [Cu; Cv];
    %disp(C)
    N = size(Cu0,1);
    iteration = 1;
    while(res > tolerance && iteration < max_iterations)
        disp(res)
        %disp(C)
        C_prev = C;
        disp(iteration)
        eq = equations(Cu, Cv, Ku, Kv, Fu, Fv1, Hu1, Hu2, Hv1, Hv2, V_mu, K_mu, K_mv, K_mfu, rq, V_mfv);
        jac = jacobian(Ku, Kv, Fu, Cu, Cv, V_mu, K_mu, K_mv, Hu1, Fv1, K_mfu, Hv1, V_mfv,rq);
        if iteration == 1
            oui = eq;
            for i = 1:numel(oui)
                %fprintf('Value at index %d: %d\n', i, oui(i));
            end
        end
        %size(C)
        %size(equations(Cu, Cv, Ku, Kv, Fu, Fv1, Hu1, Hu2, Hv1, Hv2, V_mu, K_mu, K_mv, K_mfu, rq, V_mfv))
        %size(jacobian(Ku, Kv, Fu, Cu, Cv, V_mu, K_mu, K_mv, Hu1, Fv1, K_mfu, Hv1, V_mfv,rq))
        %disp(jacobian(Ku, Kv, Fu, Cu, Cv, V_mu, K_mu, K_mv, Hu1, Fv1, K_mfu, Hv1, V_mfv,rq))
        %size(C)
        %size((jacobian(Ku, Kv, Fu, Cu, Cv, V_mu, K_mu, K_mv, Hu1, Fv1, K_mfu, Hv1, V_mfv,rq)\equations(Cu, Cv, Ku, Kv, Fu, Fv1, Hu1, Hu2, Hv1, Hv2, V_mu, K_mu, K_mv, K_mfu, rq, V_mfv)))
        C = C-(jac\eq);
        %size(C)
        Cu = C(1:N);
        Cv = C(N+1:end);
        %size(C)
        %size(Cu)
        %size(Cv)
        %disp(N)
        res = norm(C-C_prev,2);
        iteration = iteration+1;
    end
   %disp(Cu)
   %disp(Cv)
end