%% correlation for full trajectories
% h
RR = zeros(L, 1);
for t = 1:L
    RR(t) = sum(v(1:M - L).*v(t:t+M-L-1))/(M - L);
end

RR_test = zeros(L, 1);
for t = 1:L
    RR_test(t) = sum(v(1:M-t+1).*v(t:M))/(M - t+1);
end



%%
% figure;hold on
% c = xcov(v, v,'unbiased');
% cc = xcov(gradient(v, dt), v,'biased');
% plot(RR(1:1000))
% plot(RR_test(1:1000))
% plot(c(M:M+1001), 'o')
% plot(cc(M:M+1001), 'o')
% plot(h_true(tgrid_short), '*')
% legend('same len', 'diff len', 'xcor unbias', 'xcor bias', 'true')

%%
% figure;hold on
% c = xcov(v, v,'unbiased');
% cc = xcov(v, v,'biased');
% plot(RR(1:1000)/RR(1))
% plot(RR_test(1:1000)/RR_test(1))
% plot(c(M:M+1001)/c(M), 'o')
% plot(cc(M:M+1001)/cc(M), 'o')
% plot(h_true(tgrid_short), '*')
% legend('same len', 'diff len', 'xcor unbias', 'xcor bias', 'true')
% title('Fix h(0) = 1')
%%

% g
g_RR = zeros(L, 1);
for t = 1:L
    g_RR(t) = sum(v(1:M - L).*(gradient(v(t:t+M-L-1), dt) + F(v(t:t+M-L-1))))/(M-L);
end

%%
figure;hold on;
plot(gradient(RR, dt))
plot(g_RR, '.');

%%
%f
g_dif_RR = zeros(L, 1);
for t = 1:L
    g_dif_RR(t) = sum(v(1:M - L).*(F(v(t:t+M-L-1))))/(M-L);
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



%%
figure;hold on
c = xcov(v, v,'unbiased');
cc = xcov(v, v,'biased');

c_obs = xcov(v_obs, v_obs, 'unbiased');

plot(tgrid_short(1:1000), RR(1:1000), 'DisplayName', 'Same len')
plot(tgrid_short(1:1000), old_cor.RR(1:1000), 'DisplayName', 'diff len')
plot(tgrid_short(1:1000), RR_test(1:1000), 'DisplayName', 'old')
plot(tgrid_short(1:1000), c(M:M+999), 'o', 'DisplayName', 'xcor unbias')

plot(tgrid_short(1:1000), cc(M:M+999), 'o', 'DisplayName', 'xcor bias')
plot(tgrid_short(1:1000), h_true(tgrid_short(1:1000)), '*', 'DisplayName', 'true')

plot(tgrid_obs, c_obs(obs_len:end), 'O', 'DisplayName', 'obs xcor unbias')
plot(tgrid_obs_short, old_cor.RR_obs, '+', 'DisplayName', 'obs old')
xlim([0, 30])
legend()
title('Compare of RR')


%%
c_gradient = xcov(gradient(v, dt), v, 'unbiased');
gradient_c = gradient(c, dt);

figure;hold on;
plot(c_gradient(M:end), 'DisplayName', 'c grad')
plot(gradient_c(M:end), 'DisplayName', 'grad c')
plot(gradient(h_true(tgrid_short(1:1000)), dt), 'DisplayName', 'true')
xlim([0, 1000]);
legend()



% 
% 
% plot(tgrid_short(1:1000), RR(1:1000), 'DisplayName', 'Same len')
% plot(tgrid_short(1:1000), old_cor.RR(1:1000), 'DisplayName', 'diff len')
% plot(tgrid_short(1:1000), RR_test(1:1000), 'DisplayName', 'old')
% plot(tgrid_short(1:1000), c(M:M+999), 'o', 'DisplayName', 'xcor unbias')
% 
% plot(tgrid_short(1:1000), cc(M:M+999), 'o', 'DisplayName', 'xcor bias')
% plot(tgrid_short(1:1000), h_true(tgrid_short(1:1000)), '*', 'DisplayName', 'true')
% 
% plot(tgrid_obs, c_obs(obs_len:end), 'O', 'DisplayName', 'obs xcor unbias')
% plot(tgrid_obs_short, old_cor.RR_obs, '+', 'DisplayName', 'obs old')
% xlim([0, 30])
% legend()
% title('Compare of RR')



