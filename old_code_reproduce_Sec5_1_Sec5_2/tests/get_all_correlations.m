function corr_func = get_all_correlations(v, dt, F)

% M = length(v);
% corr_func.h = xcov(v, v, 'unbiased')';
% corr_func.h = corr_func.h(M:end);
% 
% corr_func.dh = gradient(corr_func.h, dt);
% corr_func.ddh = 4*del2(corr_func.h, dt);
% 
% corr_func.f = xcov(v, F(v), 'unbiased')';
% corr_func.f = corr_func.f(M:end);
% 
% corr_func.g = corr_func.dh - corr_func.f;
% corr_func.dg = corr_func.ddh - gradient(corr_func.f, dt);

M = length(v);
corr_func.h = xcov(v, v, 'unbiased')';
corr_func.h = corr_func.h(M:end);

corr_func.dh = xcov(gradient(v, dt), v, 'unbiased')';
corr_func.dh = corr_func.dh(M:end);

corr_func.ddh = xcov(4*del2(v, dt), v, 'unbiased')';
corr_func.ddh = corr_func.ddh(M:end);

corr_func.f = xcov(F(v), v, 'unbiased')';
corr_func.f = corr_func.f(M:end);

corr_func.df = xcov(gradient(F(v), dt), v, 'unbiased')';
corr_func.df = corr_func.df(M:end);

corr_func.g = -corr_func.dh + corr_func.f;
corr_func.dg = -corr_func.ddh + corr_func.df;


end