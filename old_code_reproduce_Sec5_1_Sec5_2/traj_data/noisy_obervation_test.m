obs_gap = 40;
obs_len = floor(M/obs_gap);
obs_std = 0.1;
tgrid_obs = tgrid(1:obs_gap:obs_len*obs_gap);
dt_obs = obs_gap *dt;
T_obs = obs_len*dt_obs;

v_obs = v(1:obs_gap:obs_len*obs_gap) + randn(1, obs_len)*obs_std;


v_obs_true = v(1:obs_gap:obs_len*obs_gap);
WW = v_obs - v_obs_true;


MM = length(v_obs);
T = MM*dt;
RR = zeros(MM, 1);
for t = 1:MM
    RR(t) = sum(v_obs(1:end-t+1).*v_obs(t:end))*dt/T;
end

RR_1 = zeros(MM, 1);
for t = 1:MM
    RR_1(t) = sum(v_obs_true(1:end-t+1).*v_obs_true(t:end))*dt/T;
end

RR_2 = zeros(MM, 1);
for t = 1:MM
    RR_2(t) = sum(v_obs_true(1:end-t+1).*WW(t:end))*dt/T;
end
RR_3 = zeros(MM, 1);
for t = 1:MM
    RR_3(t) = sum(WW(1:end-t+1).*v_obs_true(t:end))*dt/T;
end

RR_4 = zeros(MM, 1);
for t = 1:MM
    RR_4(t) = sum(WW(1:end-t+1).*WW(t:end))*dt/T;
end

RR_5 = RR_1 + RR_2 + RR_3 + RR_4;
figure; hold on;
plot(tgrid_obs, RR, 'LineWidth', 4)
plot(tgrid_obs, RR_1, 'red', 'LineWidth', 2)
plot(tgrid_obs, RR_2, 'o')
plot(tgrid_obs, RR_3, 'o')
plot(tgrid_obs, RR_4, '.')
% plot(tgrid_obs, RR_5)
xlim([0, 20])


%%
N = obs_len;
k_seq = 1:N;
omega = 0.0001;
figure;hold on
plot(k_seq, log10(1./(sqrt(N - k_seq))));
plot(k_seq, log10(exp(-omega*k_seq)));
plot(k_seq, log10(abs(RR - RR_1)./abs(RR_1)))
% xlim([1, 100])
