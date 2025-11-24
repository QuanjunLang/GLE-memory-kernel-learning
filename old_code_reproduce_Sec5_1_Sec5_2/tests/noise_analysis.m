close all
clear all
% rng(2)
% rng(7)

random_seed = 10;
rng(random_seed)

addpath(genpath(fileparts(which(mfilename))));
%% Generate true memory kernel along with its Fourier transform and Laplace transform
memory_para.p = 1;
memory_para.q = 3;
memory_para.plotON = 1;
[R, S, R_lap, u_seq, eta_seq, h_true] = generate_memory_kernel(memory_para);


%% Generate correlated noise and trajectories
traj_para.loadON = 1;
traj_para.wu = 80*pi;
traj_para.N = 8000;

M = 2^16;
F = @(x) 0*x;

[v, dt, tgrid, T0, correlated_noise] = generate_trajectories(R, S, F, M, traj_para, random_seed);


omega = 0.01;
rho = @(x) exp(-omega*x);        % omega
RR      = auto_correlation(v, dt);
% second derivatives of VACF
dg_true = del2(RR, dt)*4;
%% noise sequence
Len = 30;
noise_seq = 10.^linspace(-5, 1, Len);
m_h_eps_seq = zeros(Len, 1);
M_h_eps_seq = zeros(Len, 1);
m_gamma_seq = zeros(Len, 1);
M_gamma_seq = zeros(Len, 1);

theta_regression_error_seq = zeros(Len, 1);
theta_prony_error_seq = zeros(Len, 1);


h_error_seq = zeros(Len, 1);
dg_error_seq = zeros(Len, 1);
gamma_regression_error_seq = zeros(Len, 1);
gamma_prony_error_seq = zeros(Len, 1);

gamma_regression_unif_error_seq = zeros(Len, 1);
gamma_prony_unif_error_seq = zeros(Len, 1);

gamma_regression_ill_error_seq = zeros(Len, 1);
gamma_regression_ill_unif_error_seq = zeros(Len, 1);


for i = 1:length(noise_seq)
    i
    % sparse and noisy observation of v
    obs_gap = 20;
    obs_len = floor(M/obs_gap/2);
    obs_std = noise_seq(i);
    tgrid_obs = tgrid(1:obs_gap:obs_len*obs_gap);
    dt_obs = obs_gap *dt;
    T_obs = obs_len*dt_obs;
    v_obs = v(1:obs_gap:obs_len*obs_gap) + randn(1, obs_len)*obs_std;
    RR_obs  = auto_correlation(v_obs, dt_obs);
    
    
    % Prony for h
    VACF_obs_len = 40;
    Prony_I.prony_p = 20;
    Prony_I.prony_N = VACF_obs_len;
    Prony_I.obs_dx = dt_obs;
    Prony_I.obs_h_grid = RR_obs(1:VACF_obs_len);
    Prony_I.rho = rho;
    Prony_I.polycoef_method = 'LS';
    Prony_I.weight_method = 'LS';
    Prony_I.root_normalization = 1;
    Prony_I.lambda_augmentation = 1;
    
    result_prony = prony_method(Prony_I);
    
    Prony_g = result_prony.g;
    Prony_dg = result_prony.dg;
    Prony_h = result_prony.h;
    
    
    h_error_seq(i) = sqrt(sum((RR - Prony_h(tgrid')).^2.*rho(tgrid'))*dt);
    dg_error_seq(i) = sqrt(sum((dg_true - Prony_dg(tgrid')).^2.*rho(tgrid'))*dt);
    
    [m_h_eps_seq(i), M_h_eps_seq(i)] = get_m_M(omega, result_prony.h_lap);
    [m_gamma_seq(i), M_gamma_seq(i)] = get_m_M(omega, R_lap);
    
    regression_para.lb = 0;
    regression_para.rb = 30;
    regression_para.knot_num = 60;
    regression_para.deg = 3;
    regression_para.free_bdry = 1;
    regression_para.N = 4001;
    regression_para.omega = omega;
    regression_para.reg_method = 'LS';
    regression_para.h_conv_phi_dx_method = 'direct';
    h_est = @(x) real(result_prony.h(x));
    dg_est = @(x) real(result_prony.dg(x));
    g_est = @(x) real(result_prony.g(x));
    theta = get_theta_regression(regression_para, h_est, dg_est);
    
    
    result_gamma = get_gamma_from_prony_h(Prony_I, result_prony);
    theta_prony = @(x) real(result_gamma.gamma(x));
    
    
    gamma_regression_error_seq(i) = sqrt(sum((theta(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    gamma_prony_error_seq(i) = sqrt(sum((theta_prony(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    
    gamma_regression_unif_error_seq(i) = sqrt(sum((theta(tgrid') - R(tgrid')).^2)*dt);
    gamma_prony_unif_error_seq(i) = sqrt(sum((theta_prony(tgrid') - R(tgrid')).^2)*dt);
    
    theta_ill = get_theta_regression_ill_posed(regression_para, h_est, g_est);
    gamma_regression_ill_error_seq(i) = sqrt(sum((theta_ill(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    gamma_regression_ill_unif_error_seq(i) = sqrt(sum((theta_ill(tgrid') - R(tgrid')).^2)*dt);
    
end

%%
% regression_para.omega = 0;
% regression_para.reg_method = 'RKHS';
% theta = get_theta_regression(regression_para, h_est, dg_est);
% a =  sqrt(sum((theta(tgrid') - R(tgrid')).^2)*dt);
%%
figure; hold on;
plot(log10(noise_seq), log10(m_h_eps_seq), 'o', 'LineWidth', 4);
plot(log10(noise_seq), log10(M_h_eps_seq), 'o', 'LineWidth', 4);
plot(log10(noise_seq), log10(m_gamma_seq), 'LineWidth', 4);
plot(log10(noise_seq), log10(M_gamma_seq), 'LineWidth', 4);
legend('m_h', 'M_h', 'm_\gamma' ,'M_\gamma')

%%
% figure; hold on;

%%
figure;hold on;
plot(log10(noise_seq), log10(h_error_seq))
plot(log10(noise_seq), log10(dg_error_seq))
legend('h error', 'dg error')
%%
figure; hold on; grid on;
error_upper_bound_seq = 1./m_h_eps_seq.*(dg_error_seq + M_gamma_seq.*h_error_seq);
% loglog(omega_seq, upper_bound_errror, 'LineWidth', 4);
plot(log10(noise_seq), log10(error_upper_bound_seq), '-', 'Color', 'green', 'DisplayName', 'error_upper_bound', 'LineWidth', 4);
plot(log10(noise_seq), log10(1./m_h_eps_seq), '-.', 'Color', 'green', 'DisplayName', 'coercivity_constant', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_regression_error_seq), 'Color', 'black','DisplayName', 'regression error L^2(rho)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_prony_error_seq), 'Color', 'blue','DisplayName', 'Prony error L^2(rho)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_regression_unif_error_seq),'-.', 'Color', 'black','DisplayName', 'regression error L^2(Lebesgue)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_prony_unif_error_seq),  '-.', 'Color', 'blue','DisplayName', 'Prony error L^2(Lebesgue)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_regression_ill_error_seq), 'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(rho)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_regression_ill_unif_error_seq), '-.', 'Color', 'red',  'DisplayName', 'ill-posed regression error L^2(Lebesgue)', 'LineWidth', 4);
% yline(log10(a))
legend()

xlabel('log 10 omega')
ylabel('log 10 L2 error')
%%

h_norm_unif = sqrt(sum(h_true(tgrid).^2));
dg_norm_unif = sqrt(sum(dg_true.^2));
theta_norm_unif = sqrt(sum(R(tgrid).^2));

h_norm_omega = sqrt(sum(h_true(tgrid).^2.*rho(tgrid)));
dg_norm_omega = sqrt(sum(dg_true.^2.*rho(tgrid')));
theta_norm_omega = sqrt(sum(R(tgrid).^2.*rho(tgrid)));

%%
figure;hold on;
plot(log10(noise_seq), log10(h_error_seq./h_norm_omega))
plot(log10(noise_seq), log10(dg_error_seq./dg_norm_omega))
legend('h rel error', 'dg rel error')
%%
figure; hold on; grid on;
error_upper_bound_seq = 1./m_h_eps_seq.*(dg_error_seq + M_gamma_seq.*h_error_seq);
% loglog(omega_seq, upper_bound_errror, 'LineWidth', 4);
plot(log10(noise_seq), log10(error_upper_bound_seq./theta_norm_omega), '-', 'Color', 'green', 'DisplayName', 'error_upper_bound', 'LineWidth', 4);
plot(log10(noise_seq), log10(1./m_h_eps_seq./theta_norm_omega), '-.', 'Color', 'green', 'DisplayName', 'coercivity_constant', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_regression_error_seq./theta_norm_omega), 'Color', 'black','DisplayName', 'regression error L^2(rho)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_prony_error_seq./theta_norm_omega), 'Color', 'blue','DisplayName', 'Prony error L^2(rho)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_regression_unif_error_seq./theta_norm_unif),'-.', 'Color', 'black','DisplayName', 'regression error L^2(Lebesgue)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_prony_unif_error_seq./theta_norm_omega),  '-.', 'Color', 'blue','DisplayName', 'Prony error L^2(Lebesgue)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_regression_ill_error_seq./theta_norm_unif), 'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(rho)', 'LineWidth', 4);
plot(log10(noise_seq), log10(gamma_regression_ill_unif_error_seq./theta_norm_unif), '-.', 'Color', 'red',  'DisplayName', 'ill-posed regression error L^2(Lebesgue)', 'LineWidth', 4);
% yline(log10(a))
legend()

xlabel('log 10 omega')
ylabel('log 10 L2 error')


