function corr_func = get_all_correlations_new(v, dt, F, m, temporal_spacial_method, corr_method)




[L, M] = size(v);

corr_method = 'xcov';                   % xcov, manual
% temporal_spacial_method = 'M_only';     % L_only, Mixed



if M == 1
    temporal_spacial_method = 'L_only';
end

switch temporal_spacial_method
    case 'L_only'
        v_temp = v(:, 1);
        corr_func.h = compute_corr(v_temp, v_temp, corr_method);
        corr_func.dh = compute_corr(gradient(v_temp, dt), v_temp, corr_method);
        corr_func.ddh = compute_corr(4*del2(v_temp, dt), v_temp, corr_method);
        corr_func.f = compute_corr(F(v_temp), v_temp, corr_method);
        corr_func.df = compute_corr(gradient(F(v_temp), dt), v_temp, corr_method);
    case 'Mixed'
        corr_func.h     = zeros(L, 1);
        corr_func.dh    = zeros(L, 1);
        corr_func.ddh   = zeros(L, 1);
        corr_func.f     = zeros(L, 1);
        corr_func.df    = zeros(L, 1);
        for i = 1:M
            v_temp = v(:, i);
            corr_func.h     = corr_func.h + compute_corr(v_temp, v_temp, corr_method);
            corr_func.dh    = corr_func.dh+ compute_corr(gradient(v_temp, dt), v_temp, corr_method);
            corr_func.ddh   = corr_func.ddh + compute_corr(4*del2(v_temp, dt), v_temp, corr_method);
            corr_func.f     = corr_func.f + compute_corr(F(v_temp), v_temp, corr_method);
            corr_func.df    = corr_func.df + compute_corr(gradient(F(v_temp), dt), v_temp, corr_method);
        end
        corr_func.h     = corr_func.h/M;
        corr_func.dh    = corr_func.dh/M;
        corr_func.ddh   = corr_func.ddh/M;
        corr_func.f     = corr_func.f/M;
        corr_func.df    = corr_func.df/M;
    case 'M_only'
        corr_func.h     = zeros(L, 1);
        corr_func.dh    = zeros(L, 1);
        corr_func.ddh   = zeros(L, 1);
        corr_func.f     = zeros(L, 1);
        corr_func.df    = zeros(L, 1);
        for i = 1:M
            v_temp = v(:, i);
            corr_func.h     = corr_func.h + compute_corr(v_temp, v_temp, 'M_only');
            corr_func.dh    = corr_func.dh+ compute_corr(gradient(v_temp, dt), v_temp, 'M_only');
            corr_func.ddh   = corr_func.ddh + compute_corr(4*del2(v_temp, dt), v_temp, 'M_only');
            corr_func.f     = corr_func.f + compute_corr(F(v_temp), v_temp, 'M_only');
            corr_func.df    = corr_func.df + compute_corr(gradient(F(v_temp), dt), v_temp, 'M_only');
        end
        corr_func.h     = corr_func.h/M;
        corr_func.dh    = corr_func.dh/M;
        corr_func.ddh   = corr_func.ddh/M;
        corr_func.f     = corr_func.f/M;
        corr_func.df    = corr_func.df/M;

end


corr_func.g = -m*corr_func.dh + corr_func.f;
corr_func.dg = -m*corr_func.ddh + corr_func.df;

end
