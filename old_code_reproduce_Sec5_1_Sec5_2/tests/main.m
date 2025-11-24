clc
close all
clear all
% rng(2)
% rng(7)

random_seed = 15;
rng(random_seed)

addpath(genpath(fileparts(which(mfilename))));
%% Generate true memory kernel along with its Fourier transform and Laplace transform
memory_para.p = 1;
memory_para.q = 3;
memory_para.plotON = 1;
[R, S, R_lap, u_seq, eta_seq, h_true, h_true_lap] = generate_memory_kernel(memory_para);


%% Generate correlated noise and trajectories
traj_para.loadON = 0;
traj_para.wu = 80*pi;
traj_para.N = 8000;
traj_para.beta = 1;

M = 2^16;


F_type = 0;
switch F_type
    case 0
        F = @(x) 0*x;
    case 'double_well'
        E0 = 1;
        C = 0.045;
        F = @(x) E0*(C*x.^4 - x.^2);
end

figure;fplot(F)
[v, dt, tgrid, T0, correlated_noise] = generate_trajectories(R, S, F, M, traj_para, random_seed);
%% sparse and noisy observation of v
obs_gap = 40;
obs_len = floor(M/obs_gap/2);
obs_std = 1e-3;
tgrid_obs = tgrid(1:obs_gap:obs_len*obs_gap);
dt_obs = obs_gap *dt;
T_obs = obs_len*dt_obs;

v_obs = v(1:obs_gap:obs_len*obs_gap) + randn(1, obs_len)*obs_std;

figure;hold on;title('Trajectory v')
plot(tgrid, v, 'linewidth', 1);
plot(tgrid_obs, v_obs, '.', 'markersize', 4);
xline(T0, 'linewidth', 3)
xlabel('time t')
legend('True', 'obs', 'Period T0')

%% correlated noise quality check
figure;subplot(2, 2, 1);
dw = traj_para.wu/traj_para.N;
wgrid = (0:dw:traj_para.wu-dw)';
T = M*dt;
plot(wgrid, S(wgrid));xlim([0, 20]);title('Spectral function S')

subplot(2, 2, 2);hold on;
RR = zeros(M, 1);
for t = 1:M
    RR(t) = sum(correlated_noise(1:end-t+1).*correlated_noise(t:end))*dt/T;
end
plot(tgrid, R(tgrid), 'LineWidth', 4);
plot(tgrid, RR, '.', 'markersize', 15);xlim([0, 8]);
title('Autocorrelation');legend('true', 'approx');

subplot(2, 2, [3, 4]);hold on;
plot(tgrid, correlated_noise);xlabel('t');title('Trajectory');

%% Generate VACF
% True h and obs h
RR      = auto_correlation(v, dt);
RR_obs  = auto_correlation(v_obs, dt_obs);


% True g and obs g
% switch F_type
%     case 0
%         F = @(x) 0*x;
%     case 'double_well'
%         E0 = 1;
%         C = 0.045;
%         F = @(x) E0*(C*x.^4 - x.^2);
% end
% 
% 
g_traj = zeros(obs_len, 1);
for t = 1:obs_len
    g_traj(t) = -sum((gradient(v_obs(1:end-t+1), dt_obs) - F(v_obs(1:end-t+1))).*v_obs(t:end))*dt_obs/T_obs;
end
% 
% % second derivatives of VACF
dg_true = del2(RR, dt)*4;



% derivatives of VACF for later use
g_RR = gradient(RR, dt);
g_RR_obs = gradient(RR_obs, dt_obs);
%% Plot VACF
figure; hold on;
plot(tgrid, RR, 'linewidth', 2)
plot(tgrid_obs, RR_obs, 'o')
plot(tgrid, h_true(tgrid), 'red', 'linewidth', 2)
xlim([0, 20])
legend('Full VACF', 'Obs VACF', 'True h')


%%
% figure; hold on;
% plot(tgrid, g_RR, 'LineWidth', 4);
% % plot(tgrid, Prony_g(tgrid), '.', 'Linewidth', 4);
% plot(tgrid_obs, g_RR_obs, 'LineWidth', 4);
% plot(tgrid_obs, g_traj, 'LineWidth', 4);
% xlim([0, 10])
% legend('g from RR (True)', 'g from Prony (The one we should use)', 'g derived from d RR obs', 'g estimted from obs traj')

%% Measure rho

omega = 0.001;
rho = @(x) exp(-omega*x);        % omega


%% Prony method for VACF h
% Assume a finite observation of h,
% Use prony method to derive the analytical h
VACF_obs_len = 30;

Prony_I.prony_p = 6;
Prony_I.prony_N = VACF_obs_len;
Prony_I.obs_dx = dt_obs;
Prony_I.obs_h_grid = RR_obs(1:VACF_obs_len);
Prony_I.rho = rho;
Prony_I.polycoef_method = 'MP';
Prony_I.weight_method = 'LS_h0';
Prony_I.root_normalization = 1;
Prony_I.lambda_augmentation = 1;

% Prony_I_g = Prony_I;
% Prony_I_g.obs_h_grid = g_traj(1:VACF_obs_len);
% Prony_I_g.


% result = estimate_h_g_prony(Prony_I, Prony_I_g);




result_prony = prony_method(Prony_I);

Prony_g = result_prony.g;
Prony_dg = result_prony.dg;
Prony_h = result_prony.h;

% error
h_error = sqrt(sum((RR - result_prony.h(tgrid')).^2.*rho(tgrid'))*dt);
dg_error = sqrt(sum((dg_true - Prony_dg(tgrid')).^2.*rho(tgrid'))*dt);

h_error
dg_error


%% plot for Prony
% h
figure; subplot(221);hold on;
plot(tgrid, RR, 'LineWidth', 3)
plot(tgrid_obs(1:VACF_obs_len), RR_obs(1:VACF_obs_len), 'o', 'LineWidth', 3)
plot(tgrid, result_prony.h(tgrid), 'LineWidth', 3)
xlim([0, 20])
legend('True VACF', 'Obs VACF', 'Est VACF')


%g

subplot(222);hold on;
plot(tgrid, g_RR, 'LineWidth', 4);
plot(tgrid, Prony_g(tgrid), '.', 'Linewidth', 4);
plot(tgrid_obs, g_RR_obs, 'LineWidth', 4);
plot(tgrid_obs, g_traj, 'LineWidth', 4);
xlim([0, 10])
legend('g from RR (True)', 'g from Prony (The one we should use)', 'g derived from d RR obs', 'g estimted from obs traj')

%dg

subplot(223);hold on;
plot(tgrid, dg_true, 'LineWidth', 4 );
plot(tgrid, Prony_dg(tgrid), '.', 'LineWidth', 4 )
plot(tgrid_obs, gradient(g_RR_obs, dt_obs), 'LineWidth', 4 );

xlim([0, 20])
legend('dg from RR (True)', 'dg from Prony h (The one we should use)', 'dg from RR obs')


% modes
subplot(224);hold on;
scatter(real(result_prony.r), imag(result_prony.r), 130, '.'); hold on;circle(0, 0, 1)
xrange = 1.4;
xlim([-xrange, xrange]);ylim([-xrange, xrange])
title('Prony modes')


%% Estimate Gamma using Laplace transform
result_gamma = get_gamma_from_prony_h(Prony_I, result_prony);
theta_Prony = @(x) real(result_gamma.gamma(x));

%% Plot Gamma from Laplace transform
figure;
subplot(121);hold on;

gamma = R;
xgrid = tgrid_obs;
xgrid_log = 10.^(linspace(-5, 5, 500));

plot(xgrid, real(result_gamma.gamma(xgrid)), '-.', 'MarkerSize', 10, 'displayname', 'est')
plot(xgrid, gamma(xgrid), 'displayname', 'true Lap[gamma]')
legend;
ylim([-4, 4])
xlim([0, 20])
title('Memory kernel')

subplot(122);hold on;grid on;
plot(log10(xgrid_log), log10(real(result_gamma.gamma(xgrid_log))), '-.', 'MarkerSize', 10, 'displayname', 'est')
plot(log10(xgrid_log), log10(gamma(xgrid_log)), 'displayname', 'true Lap[gamma]')
legend;

ylim([-2, 4])
title('Memory kernel, log scale')
%% Using the structure with coercivity. Firstly the regression matrix,
regression_para.lb = 0;
regression_para.rb = 30;
regression_para.knot_num = 20;
regression_para.deg = 3;
regression_para.free_bdry = 1;
regression_para.N = 2001;
regression_para.omega = omega;
regression_para.reg_method = 'RKHS';
regression_para.h_conv_phi_dx_method = 'basis_derivative';
h_est = @(x) real(result_prony.h(x));
dg_est = @(x) real(result_prony.dg(x));
g_est = @(x) real(result_prony.g(x));

theta = get_theta_regression(regression_para, h_est, dg_est);
theta_ill_posed = get_theta_regression_ill_posed(regression_para, h_est, g_est);

xgrid_regression = linspace(regression_para.lb, regression_para.rb, regression_para.N)';
dx_regression = xgrid_regression(2) - xgrid_regression(1);

%% error
gamma_regression_ill_posed_error = sqrt(sum((theta_ill_posed(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)
gamma_regression_error = sqrt(sum((theta(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)
gamma_prony_error = sqrt(sum((theta_Prony(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)
%% Plot estimation result
figure;hold on;

plot(xgrid_regression, R(xgrid_regression), 'LineWidth', 4);
plot(xgrid_regression, theta(xgrid_regression), 'LineWidth', 4);

plot(xgrid_regression, real(result_gamma.gamma(xgrid_regression)), 'o', 'MarkerSize', 4)
plot(xgrid_regression, theta_ill_posed(xgrid_regression), 'LineWidth', 4);

ylim([min(theta(xgrid_regression))*1.5, max(theta(xgrid_regression))*1.5])
legend('True', 'Regression', 'Laplace transform', 'ill-posed')
title('Estimation result of memory kernel')
%% Coercivity analysis
omega = 0.001;
% Firstly, the bound for laplace transform of gamma
tau_grid_real = linspace(-20, 20, 100000);
tau_grid = omega + 1i*tau_grid_real;
z_gamma_z = abs(tau_grid.*R_lap(tau_grid));

m_gamma = min(z_gamma_z);
M_gamma = max(z_gamma_z);

% Secondly, the bound for laplace transform of h_eps
% tau_grid = omega + 1i*tau_grid_real;
z_h_eps_z = abs(tau_grid.*result_prony.h_lap(tau_grid));
h_eps_z = abs(result_prony.h_lap(tau_grid));

m_h_eps = min(z_h_eps_z.^2);
M_h_eps = max(z_h_eps_z.^2);





figure;
% subplot(131);hold on;
plot(tau_grid_real, z_gamma_z, 'LineWidth', 4)
yline(M_gamma, 'LineWidth', 4)
yline(m_gamma, 'LineWidth', 4)
legend('$|z\widehat{\gamma}(z)|$', 'm', 'M','Interpreter','latex')
xlabel('\tau')
title('The value of $|z\widehat{\gamma}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')

figure;
subplot(131);hold on;
plot(tau_grid_real, z_h_eps_z.^2, 'LineWidth', 4)
yline(M_h_eps, 'LineWidth', 4)
yline(m_h_eps, 'LineWidth', 4)
legend('$|z\widehat{h_{\varepsilon}}(z)|^2$', 'm', 'M','Interpreter','latex', 'Location', 'best')
xlabel('\tau')
title('The value of $|z\widehat{h_{\varepsilon}}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
ylim([-1, 3])


m_h_eps_z = min(h_eps_z.^2);
M_h_eps_z = max(h_eps_z.^2);

subplot(132);hold on;
plot(tau_grid_real, h_eps_z.^2, 'LineWidth', 4)
yline(M_h_eps_z, 'LineWidth', 4)
yline(m_h_eps_z, 'LineWidth', 4)

legend('$|\widehat{h_{\varepsilon}}(z)|^2$', 'm', 'M','Interpreter','latex', 'Location', 'best')
xlabel('\tau')
title('The value of $|\widehat{h_{\varepsilon}}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
ylim([-1, 3])




a1 = sum(result_prony.h(tgrid))*dt;
a2 = result_prony.h(0);
comb_coef = a1^2/(a1^2 + a2^2);
% comb_coef = 0.05;


h_eps_z_comb2 = comb_coef*z_h_eps_z.^2 + (1-comb_coef)*h_eps_z.^2;
m_h_eps_z_comb = min(h_eps_z_comb2);
M_h_eps_z_comb = max(h_eps_z_comb2);

subplot(133);hold on;
plot(tau_grid_real, h_eps_z_comb2, 'LineWidth', 4)
yline(m_h_eps_z_comb, 'LineWidth', 4)
yline(M_h_eps_z_comb, 'LineWidth', 4)
legend('$\alpha |z\widehat{h_{\varepsilon}}(z)|^2 + |\widehat{h_{\varepsilon}}(z)|^2$', 'm', 'M','Interpreter','latex', 'Location', 'best')
xlabel('\tau')
title('The value of $\alpha |z\widehat{h_{\varepsilon}}(z)|^2 + |\widehat{h_{\varepsilon}}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
ylim([-1, 3])
%%
% Loss of the regression result, using estimated h and estimated g
hgrid_regression = h_est(xgrid_regression);
theta_grid_regression = theta(xgrid_regression);
theta_grid_prony = theta_Prony(xgrid_regression);
dg_grid_regression = dg_est(xgrid_regression);
g_grid_regression = Prony_g(xgrid_regression);
rho_grid_regression = rho(xgrid_regression);

h_conv_theta = Volterra_convolution(hgrid_regression, theta_grid_regression, dx_regression);
h_conv_theta_dx = gradient(h_conv_theta, dx_regression);


h_conv_theta_prony = Volterra_convolution(hgrid_regression, theta_grid_prony, dx_regression);
h_conv_theta_prony_dx = gradient(h_conv_theta_prony, dx_regression);

h_conv_R = Volterra_convolution(hgrid_regression, R(xgrid_regression), dx_regression);
h_conv_R_dx = gradient(h_conv_R, dx_regression);



loss_est_regression = sum((h_conv_theta_dx + dg_grid_regression).^2.*rho_grid_regression)*dx_regression
loss_est_prony = sum((h_conv_theta_prony_dx + dg_grid_regression).^2.*rho_grid_regression)*dx_regression
loss_est_true = sum((h_conv_R_dx + dg_grid_regression).^2.*rho_grid_regression)*dx_regression

%% From the theorem, this is the upper bound of the error

upbd = 1/m_h_eps*dg_error + M_gamma/m_h_eps*h_error;
upbd


%%
% figure; hold on;
% plot(-g_grid_regression)
% plot(h_conv_theta_prony)
% plot(h_conv_theta)
% ylim([-10, 10])
% omega
% 


% prediction error 
% generate the trajectory v
% pred_v_regression = generate_pred_v(theta, M, dt, tgrid, F, correlated_noise, v(1));
% pred_v_Prony = generate_pred_v(theta_Prony, M, dt, tgrid, F, correlated_noise, v(1));

%%
figure; hold on;
plot(tgrid_obs, pred_v_regression(1:obs_gap:obs_len*obs_gap));
plot(tgrid_obs, pred_v_Prony(1:obs_gap:obs_len*obs_gap));
plot(tgrid_obs, v(1:obs_gap:obs_len*obs_gap))
plot(tgrid_obs, v_obs)
legend('pred traj regression', 'pred traj prony', 'True traj', 'obs traj')


traj_error_regression = sqrt(sum((pred_v_regression(1:obs_gap:obs_len*obs_gap) - v_obs).^2*dt_obs)/T_obs)
traj_error_prony = sqrt(sum((pred_v_Prony(1:obs_gap:obs_len*obs_gap) - v_obs).^2*dt_obs)/T_obs)
traj_error_true = sqrt(sum((v(1:obs_gap:obs_len*obs_gap) - v_obs).^2*dt_obs)/T_obs)

sqrt(sum((pred_v_Prony - v(1:M/2)).^2*dt)/T*2)
sqrt(sum((pred_v_regression - v(1:M/2)).^2*dt)/T*2)

sqrt(sum((v(1:M/2)).^2*dt)/T*2)
