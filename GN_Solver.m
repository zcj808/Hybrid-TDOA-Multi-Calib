function [g,norm_dk,value_f] = GN_Solver(g)
    GN_t=g.lim_t;
    if g.label=="init"
        W=g.W_init;
    elseif g.label=="final"
        W=g.W;
    end
    [ini_mic_err,ini_euler_err,ini_off_err,ini_dri_err,ini_s_err] = compute_error(g);
    [J,r] = compute_J(g);
    ini_value_f = 1/2*r'/W*r;
    tic
    while toc<GN_t/10
        dk = -(J'/W*J)\J'/W*r;
        norm_dk = norm(dk);
        xk = high2low(g);
        xk = xk+dk;
        g = low2high(xk,g);
        [J,r] = compute_J(g);
        value_f = 1/2*r'/W*r;
        if norm_dk < g.dk_p || value_f<g.f_p % use && lead gt_rec low conv. in large TDOA noise
            [mic_err,euler_err,off_err,dri_err,s_err] = compute_error(g);
            g.init_rec=[g.id,ini_value_f,ini_mic_err,ini_euler_err,ini_off_err,ini_dri_err,ini_s_err,toc];
            g.rec=[g.id,value_f,mic_err,euler_err,off_err,dri_err,s_err,toc];
            disp("convergence");
            break
        elseif value_f/ini_value_f>1e5 || isnan(value_f)
            disp("divergence");
            break
        end
    end
end