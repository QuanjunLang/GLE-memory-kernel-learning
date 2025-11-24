function corr_func = get_all_correlations_ensemble(all_v, dt, F, time_lag)
[M, N] = size(all_v);

h   = zeros(1, M);
dh  = zeros(1, M);
ddh = zeros(1, M);
f   = zeros(1, M);
df  = zeros(1, M);




% v0      = all_v(1:1+time_len, :);
all_dv = gradient(all_v', dt)';
% all_ddv = 4*del2(all_v', dt)';
all_dFv = gradient(F(all_v)', dt)';


time_lag_len = floor(time_lag/dt);

for t = 1:M
    % vt   = all_v(t, :);
    % 
    % temp = cov(v0, vt); temp = sum(v0.*vt)/(N-1);
    % h(t)    =  temp;
    % temp = my_cov(v0, gradient(vt, dt));       dh(t)   =  my_cov(v0, gradient(vt, dt))
    % temp = cov(v0, 4*del2(vt, dt));         ddh(t)  =  temp(1, 2); 
    % 
    % temp = cov(v0, F(vt));                  f(t)    =  temp(1, 2);
    % temp = cov(v0, gradient(F(vt), dt));    df(t)   =  temp(1, 2);
    
    upperbound = min(t+time_lag_len, M);
    v0      = all_v(1:1+upperbound-t, :);

    vt   = all_v(t:upperbound, :);
    dvt  = all_dv(t:upperbound, :);
    % ddvt = all_ddv(t:upperbound, :);
    dFvt = all_dFv(t:upperbound, :);

    h(t)    =  my_cov(v0, vt);
    dh(t)   =  my_cov(v0, dvt);

    % ddh(t)  =  my_cov(v0, ddvt); 

    f(t)    =  my_cov(v0, F(vt)); 
    df(t)   =  my_cov(v0, dFvt);
end

dh_gd =  gradient(h, dt);
ddh = 4*del2(h, dt);
g   = -dh + f;
dg  = -ddh + df;

corr_func.h     = h';
corr_func.dh    = dh';
corr_func.dh_gd = dh_gd';
corr_func.ddh   = ddh';
corr_func.f     = f';
corr_func.df    = df';
corr_func.g     = g';
corr_func.dg    = dg';
corr_func.g_gd  = -dh_gd' + f';
end


function c = my_cov(v1, v2)
    [a, b] = size(v1);
    c = sum(v1.*v2, 'all')/(a*b);
end