clc
close all
clear all
% rng(2)
% rng(7)

% random_seed = 111;
random_seed = 15;
rng(random_seed)

addpath(genpath(fileparts(which(mfilename))));
%% Generate true memory kernel along with its Fourier transform and Laplace transform
memory_para.p = 1;
memory_para.q = 3;
memory_para.plotON = 1;
[R, S, R_lap, u_seq, eta_seq, h_true, h_true_lap] = generate_memory_kernel(memory_para);


%% Generate correlated noise and trajectories
F_type = 'steep_double_well';
% F_type = 'single_well';
% F_type = '0';
switch F_type
    case '0'
        F = @(x) 0*x;
    case 'double_well'
        E0 = 0.02;
        F = @(x) -E0*x.*(x.^2 - 1);
    case 'single_well'
        E0 = 0.5;
        F = @(x) -E0*x;
    case 'steep_double_well'
        E0 = 0.1;
        F = @(x) -E0*x.*(x.^2 - 4);
end

traj_para.F_type = F_type;
traj_para.loadON = 1;
traj_para.wu = 80*pi;
traj_para.N = 8000;
traj_para.beta = 0.2;

M = 2^16;




figure;fplot(F)

tic
[v, dt, tgrid, T0, correlated_noise] = generate_trajectories(R, S, F, M, traj_para, random_seed);
toc
T = M*dt;
%% sparse and noisy observation of v
obs_gap = 40;
obs_len = floor(M/obs_gap/2);
obs_std = 0;
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
% figure;subplot(2, 2, 1);
% dw = traj_para.wu/traj_para.N;
% wgrid = (0:dw:traj_para.wu-dw)';
% T = M*dt;
% plot(wgrid, S(wgrid));xlim([0, 20]);title('Spectral function S')
%
% subplot(2, 2, 2);hold on;
% RR = zeros(M, 1);
% for t = 1:M
%     RR(t) = sum(correlated_noise(1:end-t+1).*correlated_noise(t:end))*dt/T;
% end
% plot(tgrid, R(tgrid), 'LineWidth', 4);
% plot(tgrid, RR, '.', 'markersize', 15);xlim([0, 8]);
% title('Autocorrelation');legend('true', 'approx');
%
% subplot(2, 2, [3, 4]);hold on;
% plot(tgrid, correlated_noise);xlabel('t');title('Trajectory');

%% Generate VACF

VACF_gen_len = 60;
tgrid_short = tgrid(1:VACF_gen_len*obs_gap);
T_short = tgrid_short(end);
% True g and obs g
RR = zeros(VACF_gen_len*obs_gap, 1);
for t = 1:VACF_gen_len*obs_gap
    RR(t) = sum(v(1:end-t+1).*v(t:end))*dt/T;
end

RR_obs = zeros(VACF_gen_len, 1);
for t = 1:VACF_gen_len
    RR_obs(t) = sum(v_obs(1:end-t+1).*v_obs(t:end))*dt_obs/T_obs;
end


% True g and obs g
g_RR = zeros(VACF_gen_len*obs_gap, 1);
for t = 1:VACF_gen_len*obs_gap
    g_RR(t) = -sum((gradient(v(1:end-t+1), dt) - F(v(1:end-t+1))).*v(t:end))*dt/T;
end

g_RR_obs = zeros(VACF_gen_len, 1);
for t = 1:VACF_gen_len
    g_RR_obs(t) = -sum((gradient(v_obs(1:end-t+1), dt_obs) - F(v_obs(1:end-t+1))).*v_obs(t:end))*dt_obs/T_obs;
end


% True dif_g and obs dif_g
g_dif_RR = zeros(VACF_gen_len*obs_gap, 1);
for t = 1:VACF_gen_len*obs_gap
    g_dif_RR(t) = sum(F(v(1:end-t+1)).*v(t:end))*dt/T;
end

g_dif_RR_obs = zeros(VACF_gen_len, 1);
for t = 1:VACF_gen_len
    g_dif_RR_obs(t) = sum(F(v_obs(1:end-t+1)).*v_obs(t:end))*dt_obs/T_obs;
end
% second derivatives of VACF
dg_true = gradient(g_RR, dt);
% True h and obs h
% RR      = auto_correlation(v, dt);
% RR_obs  = auto_correlation(v_obs, dt_obs);


%%
figure;
subplot(131);hold on
plot(tgrid_short(1:VACF_gen_len*obs_gap), g_RR(1:VACF_gen_len*obs_gap), 'LineWidth', 3)
plot(tgrid_obs(1:VACF_gen_len), g_RR_obs(1:VACF_gen_len), 'o', 'LineWidth', 3)
% plot(tgrid_obs(1:VACF_gen_len), g_RR_obs(1:VACF_gen_len) + g_dif_RR_obs(1:VACF_gen_len), 'o', 'LineWidth', 3)

legend('true','obs')
xlim([0, 10])

subplot(132); hold on
plot(tgrid_short(1:VACF_gen_len*obs_gap), RR(1:VACF_gen_len*obs_gap), 'LineWidth', 3)
plot(tgrid_obs(1:VACF_gen_len), RR_obs(1:VACF_gen_len), 'o', 'LineWidth', 3)
legend('true','obs')
xlim([0, 10])

subplot(133); hold on
plot(tgrid_short(1:VACF_gen_len*obs_gap), g_dif_RR(1:VACF_gen_len*obs_gap), 'LineWidth', 3)
plot(tgrid_obs(1:VACF_gen_len), g_dif_RR_obs(1:VACF_gen_len), 'o', 'LineWidth', 3)
legend('true','obs')
xlim([0, 10])

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
Prony_I.weight_method = 'LS_h0_new';
Prony_I.root_normalization = 1;
Prony_I.lambda_augmentation = 1;
Prony_I.drop_0 = 0;

% For the forcing term only
% Combine the forcing term and the derivative of h

Prony_I_g = Prony_I;
Prony_I_g.weight_method = 'LS';
Prony_I_g.obs_h_grid = g_dif_RR_obs(1:VACF_obs_len);
Prony_I_g.prony_p = Prony_I.prony_p+2;
Prony_I_g.drop_0 = 0;

[result_h, result_g, result_f] = estimate_h_g_prony(Prony_I, Prony_I_g, F_type);




% result_prony = prony_method(Prony_I);
%
Prony_g = result_g.g;
Prony_dg = result_g.dg;
Prony_h = result_h.h;

% error
h_error = sqrt(sum((RR - Prony_h(tgrid_short')).^2.*rho(tgrid_short'))*dt);
dg_error = sqrt(sum((dg_true - Prony_dg(tgrid_short')).^2.*rho(tgrid_short'))*dt);

h_error
dg_error


%% plot for Prony
% h
figure; subplot(221);hold on;
plot(tgrid_short, RR, 'LineWidth', 3)
plot(tgrid_obs(1:VACF_obs_len), RR_obs(1:VACF_obs_len), 'o', 'LineWidth', 3)
plot(tgrid_short, Prony_h(tgrid_short), 'LineWidth', 3)
xlim([0, 20])
legend('True VACF', 'Obs VACF', 'Est VACF')


%g

subplot(222);hold on;
plot(tgrid_short, g_RR, 'LineWidth', 4);
plot(tgrid_short, Prony_g(tgrid_short), '.', 'Linewidth', 4);
plot(tgrid_obs(1:VACF_gen_len), g_RR_obs, 'LineWidth', 4);
% plot(tgrid_obs(1:VACF_gen_len), result., 'LineWidth', 4);
xlim([0, 10])
legend('g from RR (True)', 'g from Prony (The one we should use)', 'g derived from d RR obs', 'g estimted from obs traj')

%dg

subplot(223);hold on;
plot(tgrid_short, dg_true, 'LineWidth', 4 );
plot(tgrid_short, Prony_dg(tgrid_short), '.', 'LineWidth', 4 )
plot(tgrid_obs(1:VACF_gen_len), gradient(g_RR_obs, dt_obs), 'LineWidth', 4 );

xlim([0, 20])
legend('dg from RR (True)', 'dg from Prony h (The one we should use)', 'dg from RR obs')


% modes
subplot(224);hold on;
scatter(real(result_h.r), imag(result_h.r), 130, '.'); hold on;circle(0, 0, 1)
xrange = 1.4;
xlim([-xrange, xrange]);ylim([-xrange, xrange])
title('Prony modes')


if F_type ~= '0'
    figure;
    
    subplot(221);hold on;
    plot(tgrid_short, g_dif_RR, 'LineWidth', 4 );
    plot(tgrid_short, result_f.h(tgrid_short), 'o', 'LineWidth', 2 )
    plot(tgrid_obs(1:VACF_gen_len), g_dif_RR_obs, 'LineWidth', 4 );
    
    xlim([0, 20])
    legend('dg from RR (True)', 'dg from Prony h (The one we should use)', 'dg from RR obs')
    
    subplot(222);hold on;
    plot(tgrid_short, gradient(g_dif_RR, dt), 'LineWidth', 4 );
    plot(tgrid_short, result_f.g(tgrid_short), 'o', 'LineWidth', 2 )
    % plot(tgrid_obs(1:VACF_gen_len), g_dif_RR_obs, 'LineWidth', 4 );
    
    xlim([0, 20])
    legend('dg from RR (True)', 'dg from Prony h (The one we should use)', 'dg from RR obs')
    
    subplot(223);hold on;
    plot(tgrid_short, gradient(gradient(RR, dt), dt), 'LineWidth', 4 );
    plot(tgrid_short, result_h.dg(tgrid_short), 'o', 'LineWidth', 2 )
    % plot(tgrid_obs(1:VACF_gen_len), g_dif_RR_obs, 'LineWidth', 4 );
    
    xlim([0, 20])
    legend('dg from RR (True)', 'dg from Prony h (The one we should use)', 'dg from RR obs')
    
end
%% Estimate Gamma using Laplace transform
switch F_type
    case '0'
        result_gamma = get_gamma_from_prony_h(Prony_I, result_h);
        theta_Prony = @(x) real(result_gamma.gamma(x));
    otherwise
        result_gamma = get_gamma_from_prony_h_esternal_force(Prony_I, Prony_I_g, result_h, result_f);
        theta_Prony = @(x) real(result_gamma.gamma(x));
end

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
regression_para.knot_num = 50;
regression_para.deg = 3;
regression_para.free_bdry = 1;
regression_para.N = 2001;
regression_para.omega = omega;
regression_para.reg_method = 'LS';
regression_para.h_conv_phi_dx_method = 'direct';
regression_para.basis_type = 'prony';
regression_para.lam = result_gamma.lam_seq;

h_est = @(x) real(Prony_h(x));
dg_est = @(x) real(Prony_dg(x));
g_est = @(x) real(Prony_g(x));

theta = get_theta_regression(regression_para, h_est, dg_est);
theta_ill_posed = get_theta_regression_ill_posed(regression_para, h_est, g_est);

xgrid_regression = linspace(regression_para.lb, regression_para.rb, regression_para.N)';
dx_regression = xgrid_regression(2) - xgrid_regression(1);

%% error
gamma_regression_ill_posed_error = sqrt(sum((theta_ill_posed(tgrid_short') - R(tgrid_short')).^2.*rho(tgrid_short'))*dt)
gamma_regression_error = sqrt(sum((theta(tgrid_short') - R(tgrid_short')).^2.*rho(tgrid_short'))*dt)
gamma_prony_error = sqrt(sum((theta_Prony(tgrid_short') - R(tgrid_short')).^2.*rho(tgrid_short'))*dt)
%% Plot estimation result
figure;hold on;

plot(xgrid_regression, R(xgrid_regression), 'LineWidth', 4);
plot(xgrid_regression, theta(xgrid_regression), 'LineWidth', 4);

plot(xgrid_regression, real(result_gamma.gamma(xgrid_regression)), 'o', 'MarkerSize', 4)
plot(xgrid_regression, theta_ill_posed(xgrid_regression), 'LineWidth', 4);

ylim([min(R(xgrid_regression))*1.5, max(R(xgrid_regression))*1.5])
legend('True', 'Regression', 'Laplace transform', 'ill-posed')
title('Estimation result of memory kernel')
%% Coercivity analysis
omega = 0.01;
% Firstly, the bound for laplace transform of gamma
tau_grid_real = linspace(-20, 20, 100000);
tau_grid = omega + 1i*tau_grid_real;
z_gamma_z = abs(tau_grid.*R_lap(tau_grid));
gamma_eps_z = abs(R_lap(tau_grid));

m_gamma = min(z_gamma_z);
M_gamma = max(z_gamma_z);

% Secondly, the bound for laplace transform of h_eps
% tau_grid = omega + 1i*tau_grid_real;
z_h_eps_z = abs(tau_grid.*result_h.h_lap(tau_grid));
h_eps_z = abs(result_h.h_lap(tau_grid));

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

y_up = 1;
figure;
subplot(131);hold on;
plot(tau_grid_real, z_h_eps_z.^2, 'LineWidth', 4)
yline(M_h_eps, 'LineWidth', 4)
yline(m_h_eps, 'LineWidth', 4)
legend('$|z\widehat{h_{\varepsilon}}(z)|^2$', 'm', 'M','Interpreter','latex', 'Location', 'best')
xlabel('\tau')
title('The value of $|z\widehat{h_{\varepsilon}}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
% ylim([-1, y_up])


m_h_eps_z = min(h_eps_z.^2);
M_h_eps_z = max(h_eps_z.^2);

subplot(132);hold on;
plot(tau_grid_real, h_eps_z.^2, 'LineWidth', 4)
yline(M_h_eps_z, 'LineWidth', 4)
yline(m_h_eps_z, 'LineWidth', 4)

legend('$|\widehat{h_{\varepsilon}}(z)|^2$', 'm', 'M','Interpreter','latex', 'Location', 'best')
xlabel('\tau')
title('The value of $|\widehat{h_{\varepsilon}}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
% ylim([-1, y_up])




a1 = sum(result_h.h(tgrid_short))*dt;
a2 = result_h.h(0);
comb_coef = a1^2/(a1^2 + a2^2);
% comb_coef = 1;
% comb_coef = 0.19;
% comb_coef = 0.2;



h_eps_z_comb2 = comb_coef*z_h_eps_z.^2 + (1-comb_coef)*h_eps_z.^2;
m_h_eps_z_comb = min(h_eps_z_comb2);
M_h_eps_z_comb = max(h_eps_z_comb2);


gamma_eps_z_comb2 = comb_coef*z_gamma_z.^2 + (1-comb_coef)*gamma_eps_z.^2;
m_gamma_eps_z_comb = min(gamma_eps_z_comb2);
M_gamma_eps_z_comb = max(gamma_eps_z_comb2);


subplot(133);hold on;
plot(tau_grid_real, h_eps_z_comb2, 'LineWidth', 4)
yline(m_h_eps_z_comb, 'LineWidth', 4)
yline(M_h_eps_z_comb, 'LineWidth', 4)
legend('$\alpha |z\widehat{h_{\varepsilon}}(z)|^2 + |\widehat{h_{\varepsilon}}(z)|^2$', 'm', 'M','Interpreter','latex', 'Location', 'best')
xlabel('\tau')
title('The value of $\alpha |z\widehat{h_{\varepsilon}}(z)|^2 + |\widehat{h_{\varepsilon}}(z)|^2$ when $z = \omega + i \tau$','Interpreter','latex')
% ylim([-1, y_up])
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

m_h_eps_z_comb
upbd_iprv = sqrt(2*(1/m_h_eps_z_comb*dg_error.^2 + M_gamma_eps_z_comb/m_h_eps_z_comb*h_error.^2));
upbd_iprv
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
% pred_v_regression = generate_pred_v(theta, M, dt, tgrid_short, F, correlated_noise, v(1));
% pred_v_Prony = generate_pred_v(theta_Prony, M, dt, tgrid_short, F, correlated_noise, v(1));

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
