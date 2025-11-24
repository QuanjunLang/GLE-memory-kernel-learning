%% correlation for full trajectories
% h


RR = zeros(L, 1);
for t = 1:L
    RR(t) = sum(v(1:M-t+1).*v(t:M))/(M - t+1);
end

%%

% g
g_RR = zeros(L, 1);
for t = 1:L
    g_RR(t) = sum(v(1:M-t+1).*(gradient(v(t:M), dt) - F(v(t:M))))/(M - t+1);
end


%%
figure;hold on;
plot(gradient(RR, dt))
plot(g_RR);


figure;plot(c(M+1:M+101));hold on;plot(RR(1:100))
%%

%f
g_dif_RR = zeros(L, 1);
for t = 1:L
    g_dif_RR(t) = sum(v(1:M - L).*(F(v(t:t+M-L-1))))/(M - t+1);
end

dg_true = gradient(g_RR, dt);

%% correlation for observed trajectories
RR_obs = zeros(L_obs, 1);
for t = 1:L_obs
    RR_obs(t) = sum(v_obs(1:M_obs - L_obs).*v_obs(t:t+M_obs-L_obs-1))/(M_obs - L_obs);
end

%g
g_RR_obs = zeros(L_obs, 1);
for t = 1:L_obs
    g_RR_obs(t) = sum(v_obs(1:M_obs - L_obs).*(-gradient(v_obs(t:t+M_obs-L_obs-1), dt_obs) + F(v_obs(t:t+M_obs-L_obs-1))))/(M_obs-L_obs);
end

%f
g_dif_RR_obs = zeros(L_obs, 1);
for t = 1:L_obs
    g_dif_RR_obs(t) = sum(v_obs(1:M_obs - L_obs).*(F(v_obs(t:t+M_obs-L_obs-1))))/(M_obs-L_obs);
end




correlations.h_full
correlations.dh_full
correlations.f_full
correlations.g_full
correlations.g_full_using_dh
