function corr_func = get_all_correlations_ensemble_mixed(all_v, dt, F)

[M, N] = size(all_v);

all_h = zeros(N, M);
all_dh = zeros(N, M);
% all_ddh = zeros(N, M);
all_f = zeros(N, M);
% all_df = zeros(N, M);









for i = 1:N
    v = all_v(:, i);
    h = xcorr(v, v, 'unbiased')';
    all_h(i, :) = h(M:end);

    f = xcorr(F(v), v, 'unbiased')';
    all_f(i, :) = f(M:end);

    dh = xcorr(gradient(v, dt), v, 'unbiased')';
    all_dh(i, :) = dh(M:end);

    % corr_func.dh = xcov(gradient(v, dt), v, 'unbiased')';
    % corr_func.dh = corr_func.dh(M:end);
    % 
    % corr_func.ddh = xcov(4*del2(v, dt), v, 'unbiased')';
    % corr_func.ddh = corr_func.ddh(M:end);

    % corr_func.f = xcov(F(v), v, 'unbiased')';
    % corr_func.f = corr_func.f(M:end);
    % 
    % corr_func.df = xcov(gradient(F(v), dt), v, 'unbiased')';
    % corr_func.df = corr_func.df(M:end);
    % 
    % corr_func.g = -corr_func.dh + corr_func.f;
    % corr_func.dg = -corr_func.ddh + corr_func.df;

end

corr_func.h = mean(all_h)';
corr_func.dh = gradient(corr_func.h, dt);
corr_func.dh_gd = mean(all_dh)';
corr_func.ddh = 4*del2(corr_func.h, dt);
corr_func.f = mean(all_f)';
corr_func.df = gradient(corr_func.f, dt);
corr_func.g = -corr_func.dh + corr_func.f;
corr_func.g_gd = -corr_func.dh_gd + corr_func.f;
corr_func.dg = -corr_func.ddh + corr_func.df;




end