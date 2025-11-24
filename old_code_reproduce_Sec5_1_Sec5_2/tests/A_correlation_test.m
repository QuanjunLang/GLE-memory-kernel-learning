%% Generate VACF
L_obs = 60;
L = L_obs*obs_gap;

tgrid_short = tgrid(1:L_obs*obs_gap);
tgrid_obs_short = tgrid_obs(1:L_obs);

%% Correlation functions of true trajectories
RR = zeros(L, 1);

for t = 1:L
    RR(t) = sum(v(1:M - L).*v(t:t+M-L-1))/(M - L);
end

dv_RR_1 = zeros(L, 1);
dv_RR_2 = zeros(L, 1);

for t = 1:L
    v_t1 = v(t+1:t+M-L);
    v_t0 = v(t:t+M-L-1);
    dv = (v_t1 - v_t0)/dt;
    dv_RR_1(t) = sum(v(1:M - L).*dv)/(M-L);
    dv_RR_2(t) = sum(v(1:M - L).*gradient(v(t:t+M-L-1), dt))/(M-L);
end


% True g and obs g
g_RR = zeros(L, 1);
for t = 1:L
    g_RR(t) = sum(v(1:M - L).*(gradient(v(t:t+M-L-1), dt) - F(v(t:t+M-L-1))))/(M-L);
end


% True dif_g and obs dif_g
g_dif_RR = zeros(L, 1);
for t = 1:L
    g_dif_RR(t) = sum(v(1:M - L).*(-F(v(t:t+M-L-1))))/(M-L);
end


g_dif_RR_wrong = zeros(VACF_gen_len*obs_gap, 1);
for t = 1:VACF_gen_len*obs_gap
    g_dif_RR_wrong(t) = sum(F(v(1:end-t+1)).*v(t:end))*dt/T;
end

% dv_RR = zeros(VACF_gen_len*obs_gap, 1);
% for t = 1:VACF_gen_len*obs_gap
%     dv_RR(t) = sum((v(2:end-t)-v(1:end-1+t)/dt).*v(t:end))*dt/T;
% end
% 
% 
% RR_obs = zeros(L_obs, 1);
% for t = 1:L_obs
%     RR_obs(t) = sum(v_obs(1:end-t+1).*v_obs(t:end))*dt_obs/T_obs;
% end
% 
% 
% g_RR_obs = zeros(L_obs, 1);
% for t = 1:L_obs
%     g_RR_obs(t) = -sum((gradient(v_obs(1:end-t+1), dt_obs) - F(v_obs(1:end-t+1))).*v_obs(t:end))*dt_obs/T_obs;
% end
% 
% 
% 
% 
g_dif_RR_obs = zeros(L_obs, 1);
for t = 1:L_obs
    g_dif_RR_obs(t) = sum(F(v_obs(1:end-t+1)).*v_obs(t:end))*dt_obs/T_obs;
end
% % second derivatives of VACF
% dg_true = gradient(g_RR, dt);

%%
RR = zeros(L_obs*obs_gap, 1);
L = L_obs*obs_gap;

for t = 1:L
    RR(t) = sum(v(1:M - L).*v(t:t+M-L-1))/(M - L);
end

dv_RR_1 = zeros(L, 1);
dv_RR_2 = zeros(L, 1);

for t = 1:L
    v_t1 = v(t+1:t+M-L);
    v_t0 = v(t:t+M-L-1);
    dv = (v_t1 - v_t0)/dt;
    dv_RR_1(t) = sum(v(1:M - L).*dv)/(M-L);
    dv_RR_2(t) = sum(v(1:M - L).*gradient(v(t:t+M-L-1), dt))/(M-L);
end

%%
figure;
subplot(221);hold on;
plot(RR);

subplot(222);hold on;
plot(gradient(RR, dt));
plot(dv_RR_1, '.')
plot(dv_RR_2)

subplot(223);hold on;
plot(g_RR);

subplot(224);

%%
figure; hold on;
plot(tgrid_short, g_RR - dv_RR_2);
plot(tgrid_short, g_dif_RR)
plot(tgrid_short, -g_dif_RR_wrong)
plot(tgrid_obs_short, -g_dif_RR_obs)

legend('true', 'another true', 'wrong full','wrong obs')