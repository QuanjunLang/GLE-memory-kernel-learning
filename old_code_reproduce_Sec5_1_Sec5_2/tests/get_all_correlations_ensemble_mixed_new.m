function corr_func = get_all_correlations_ensemble_mixed_new(all_v, dt, F)

[M, N] = size(all_v);

all_dv = gradient(all_v', dt)';
all_ddv = 4*del2(all_v', dt)';
all_dFv = gradient(F(all_v)', dt)';


h = zeros(1, M);
dh = zeros(1, M);
ddh = zeros(1, M);
f = zeros(1, M);
df = zeros(1, M);
g = zeros(1, M);
dg = zeros(1, M);

tic
parfor i = 1:M
    i
    v0_ensemble = all_v(1:M-i+1, :);
    v1_ensemble = all_v(i:M, :);
    Fv1_ensemble = F(v1_ensemble);
    % dv1_ensemble = all_dv(i:M, :);
    % ddv1_ensemble = all_ddv(i:M, :);
    
    % dFv1_ensemble = all_dFv(1:i:M, :);

    h(i) = my_cov(v0_ensemble, v1_ensemble);
    % dh(i)= my_cov(v0_ensemble, dv1_ensemble);
    % ddh(i)= my_cov(v0_ensemble, ddv1_ensemble);
    f(i) = my_cov(v0_ensemble, Fv1_ensemble);

    % f(i) = my_cov_ensemble(v0_ensemble, Fv1_ensemble);

    % df(i) = my_cov(v0_ensemble, dFv1_ensemble);

end

toc
corr_func.h = h;
corr_func.dh = gradient(corr_func.h, dt);
corr_func.ddh = 4*del2(corr_func.h, dt);
corr_func.f = f;

corr_func.df = gradient(corr_func.f, dt);


corr_func.g = -corr_func.dh + corr_func.f;
corr_func.dg = -corr_func.ddh + corr_func.df;
end


function result = my_cov(v0, v1)
[M_small, N] = size(v0);
result = sum(v1.*v0, 'all')/(M_small*N);
end

% function c = my_cov_ensemble(v1, v2)
%     c = sum(v1.*v2)/(length(v1)-1);
% end