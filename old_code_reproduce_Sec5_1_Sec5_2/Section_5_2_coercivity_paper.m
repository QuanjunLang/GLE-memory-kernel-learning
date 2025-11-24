%%
clc
close all
clear all

temp = matlab.desktop.editor.getActiveFilename;
parent_dir_temp = strsplit(temp, '/');
parent_dir = strjoin(parent_dir_temp(1:end-3), filesep);
fig_dir = [parent_dir, '/note/figure/'];
random_seed = 15;
rng(random_seed)
addpath(genpath(fileparts(which(mfilename))));


memory_para.p = 1;
memory_para.q = 3;
memory_para.plotON = 1;
[R, S, R_lap, u_seq, eta_seq, h_true, h_true_lap, dh_true, ddh_true] = generate_memory_kernel(memory_para);

% Generate correlated noise and trajectories
F_type = '0';
F = @(x) 0*x;

traj_para.F_type = F_type;
traj_para.loadON = 1;
traj_para.wu = 80*pi;
traj_para.N = 8000;
traj_para.beta = 1;

M = 2^16;

tic;[v, dt, tgrid, T0, correlated_noise] = generate_trajectories(R, S, F, M, traj_para, random_seed);toc;
T = M*dt;






% sparse and noisy observation of v%%%%%
obs_gap = 30;
obs_len = floor((M/obs_gap)/2);
obs_std = 1e-1;
VACF_obs_len = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Measure rho

omega = 0.05;
rho = @(x) exp(-2*omega*x);

tgrid_obs = tgrid(1:obs_gap:obs_len*obs_gap);
dt_obs = obs_gap *dt;
T_obs = obs_len*dt_obs;
v_obs = v(1:obs_gap:obs_len*obs_gap) + randn(1, obs_len)*obs_std;
M_obs = obs_len;

% Generate VACF
corr_func = get_all_correlations(v, dt, F);
corr_func_obs = get_all_correlations(v_obs, dt_obs, F);
% Prony method for VACF h
Prony_I.prony_p = 10;
Prony_I.prony_N = VACF_obs_len;
Prony_I.obs_dx = dt_obs;
Prony_I.obs_h_grid = corr_func_obs.h(1:VACF_obs_len);
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
Prony_I_g.obs_h_grid = corr_func_obs.f(1:VACF_obs_len);
Prony_I_g.prony_p = Prony_I.prony_p+2;
Prony_I_g.drop_0 = 0;
[result_h, result_g, result_f] = estimate_h_g_prony(Prony_I, Prony_I_g, F_type);
% comb_coef is alpha_1 in the paper
a2 = real(-sum(result_h.w./result_h.lam));
a1 = result_h.h(0);
comb_coef = a1^2/(a1^2 + a2^2);



%% Omega sequence
Len = 20;
omega_seq = 10.^linspace(-3, 2, Len);

m_h_eps_seq = zeros(Len, 1);
M_h_eps_seq = zeros(Len, 1);
m_gamma_seq = zeros(Len, 1);
M_gamma_seq = zeros(Len, 1);

gamma_regression_com_error_seq = zeros(Len, 1);
gamma_regression_dev_error_seq = zeros(Len, 1);
gamma_regression_ill_error_seq = zeros(Len, 1);
gamma_prony_error_seq = zeros(Len, 1);

h_error_seq_square = zeros(Len, 1);
g_error_seq_square = zeros(Len, 1);

h_norm_seq = zeros(Len, 1);
g_norm_seq = zeros(Len, 1);
gamma_norm_seq = zeros(Len, 1);


[result_h, result_g, result_f] = estimate_h_g_prony(Prony_I, Prony_I_g, F_type);
Prony_g = result_g.g;
Prony_dg = result_g.dg;
Prony_h = result_h.h;
h_est = @(x) real(Prony_h(x));
dg_est = @(x) real(Prony_dg(x));
g_est = @(x) real(Prony_g(x));

result_gamma = get_gamma_from_prony_h(Prony_I, result_h);
theta_Prony = @(x) real(result_gamma.gamma(x));
for i = 1:length(omega_seq)
    i
    omega = omega_seq(i);
    rho = @(x) exp(-2*omega*x);        % omega

    h_error_seq_square(i) = sum((h_true(tgrid') - Prony_h(tgrid')).^2.*rho(tgrid'))*dt;
    g_error_seq_square(i) = comb_coef*sum((-dh_true(tgrid') - Prony_g(tgrid')).^2.*rho(tgrid'))*dt + (1-comb_coef)*sum((-ddh_true(tgrid') - Prony_dg(tgrid')).^2.*rho(tgrid'))*dt;

    % bounds for gamma and h lap
    tau_grid_real = 10.^linspace(-10, 10, 10000);
    tau_grid = omega + 1i*tau_grid_real;

    z_gamma_z = abs(tau_grid.*R_lap(tau_grid));             gamma_z = abs(R_lap(tau_grid));
    gamma_z_comb = comb_coef*gamma_z.^2 + (1-comb_coef)*z_gamma_z.^2;
    m_gamma = min(gamma_z_comb);                            M_gamma = max(gamma_z_comb);

    z_h_eps_z = abs(tau_grid.*result_h.h_lap(tau_grid));    h_eps_z = abs(result_h.h_lap(tau_grid));
    h_eps_z_comb = comb_coef*h_eps_z.^2 + (1-comb_coef)*z_h_eps_z.^2;
    m_h_eps_z_comb = min(h_eps_z_comb);                     M_h_eps_z_comb = max(h_eps_z_comb);
    m_2 = min(z_h_eps_z.^2);

    m_h_eps_seq(i) = m_h_eps_z_comb;
    M_h_eps_seq(i) = M_h_eps_z_comb;
    m_2_seq(i) = m_2;
    m_gamma_seq(i) = m_gamma;
    M_gamma_seq(i) = M_gamma;

    regression_para.lb = 0;
    regression_para.rb = 30;
    regression_para.knot_num = 30;
    regression_para.deg = 3;
    regression_para.free_bdry = 1;
    regression_para.N = 2001;
    regression_para.omega = omega;
    regression_para.reg_method = 'RKHS';
    regression_para.h_conv_phi_dx_method = 'direct';
    % regression_para.basis_type = 'prony';
    regression_para.basis_type = 'spline';
    regression_para.lam = result_gamma.lam_seq;
    regression_para.comb_coef = comb_coef;


    theta_derivative = get_theta_regression(regression_para, h_est, dg_est);
    theta_ill_posed = get_theta_regression_ill_posed(regression_para, h_est, g_est);
    [theta, result_details_regression] = get_theta_regression_combined(regression_para, h_est, dg_est, g_est);

    gamma_regression_ill_posed_error = sqrt(sum((theta_ill_posed(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    gamma_regression_derivative_error = sqrt(sum((theta_derivative(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    gamma_prony_error = sqrt(sum((theta_Prony(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    gamma_regression_combined = sqrt(sum((theta(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);

    gamma_regression_com_error_seq(i) = gamma_regression_combined;
    gamma_regression_dev_error_seq(i) = gamma_regression_derivative_error;
    gamma_regression_ill_error_seq(i) = gamma_regression_ill_posed_error;
    gamma_prony_error_seq(i) = gamma_prony_error;

    h_norm_seq(i) = sqrt(sum(h_true(tgrid).^2.*rho(tgrid)));
    g_norm_seq(i) = comb_coef*sum((corr_func.g).^2.*rho(tgrid'))*dt + (1-comb_coef)*sum((corr_func.dg).^2.*rho(tgrid'))*dt;
    gamma_norm_seq(i) = sqrt(sum(R(tgrid).^2.*rho(tgrid)));

end

error_upper_bound_seq = sqrt(2./m_h_eps_seq.*(g_error_seq_square + M_gamma_seq.*h_error_seq_square));
%% color
blue = [0    0.4470    0.7410];
red =     [0.8500    0.3250    0.0980];
yellow =     [0.9290    0.6940    0.1250];
purple =    [ 0.4940    0.1840    0.5560];
green =   [  0.4660    0.6740    0.1880];

%%
fig = figure;
t = tiledlayout(1, 3,"TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';





% First graph: L2 error change with omega
nexttile;grid on;hold on;
% 
% plot(log10(omega_seq), log10(gamma_regression_com_error_seq), '-', 'Color', [ 0.8500    0.3250    0.0980], 'DisplayName', '$\|\theta - \gamma\|_{L^2(\rho)}$', 'LineWidth', 6);
% plot(log10(omega_seq), log10(gamma_prony_error_seq), ':',  'Color',  [0.9290    0.6940    0.1250], 'DisplayName', '$\|\theta_L - \gamma\|_{L^2(\rho)}$', 'LineWidth', 4, 'MarkerSize',5);
% plot(log10(omega_seq), log10(gamma_regression_ill_error_seq), '-.', 'Color', [0.4660    0.6740    0.1880], 'DisplayName', '$\|\theta_1 - \gamma\|_{L^2(\rho)}$', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_regression_dev_error_seq), '--', 'Color', [0.4940    0.1840    0.5560], 'DisplayName', '$\|\theta_2 - \gamma\|_{L^2(\rho)}$', 'LineWidth', 4);
plot(log10(omega_seq), log10(error_upper_bound_seq), '-', 'Color', [0    0.4470    0.7410],'DisplayName', ['Theoretical' newline 'upper bound'], 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_com_error_seq), '-', 'Color', [ 0.8500    0.3250    0.0980], 'DisplayName', '$\theta$', 'LineWidth', 6);
plot(log10(omega_seq), log10(gamma_regression_ill_error_seq), '-.', 'Color', [0.4660    0.6740    0.1880], 'DisplayName', '$\theta_1$', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_dev_error_seq), '--', 'Color', [0.4940    0.1840    0.5560], 'DisplayName', '$\theta_2$', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_prony_error_seq), ':',  'Color',  [0.9290    0.6940    0.1250], 'DisplayName', '$\theta_L$', 'LineWidth', 4, 'MarkerSize',5);

xlabel('$log_{10}(\omega)$','Interpreter','latex')
ylabel('$log_{10} L^2(\rho)$ error','Interpreter','latex')
xlim([min(log10(omega_seq)), max(log10(omega_seq))])
title('Change of the $L^2(\rho)$ errors with $\omega$','Interpreter','latex')
legend('Interpreter','latex')






% Second graph: Relative L2 error change with omega
nexttile;grid on;hold on;
plot(log10(omega_seq), log10(gamma_norm_seq), '-', 'Color', blue, 'DisplayName', '$\|\gamma\|_{L^2(\rho)}$', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_com_error_seq./gamma_norm_seq), '-', 'Color', red, 'DisplayName', '$\theta$', 'LineWidth', 6);
plot(log10(omega_seq), log10(gamma_regression_ill_error_seq./gamma_norm_seq), '-.', 'Color', [0.4660    0.6740    0.1880], 'DisplayName', '$\theta_1$', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_regression_dev_error_seq./gamma_norm_seq), '--', 'Color', [0.4940    0.1840    0.5560], 'DisplayName', '$\theta_2$', 'LineWidth', 4);
plot(log10(omega_seq), log10(gamma_prony_error_seq./gamma_norm_seq), ':',  'Color',  [0.9290    0.6940    0.1250], 'DisplayName', '$\theta_L$', 'LineWidth', 4, 'MarkerSize',5);

xlabel('$log_{10}(\omega)$','Interpreter','latex')
ylabel('$log_{10}$ relative $L^2(\rho)$ error','Interpreter','latex')
xlim([min(log10(omega_seq)), max(log10(omega_seq))])
title('Change of the relative $L^2(\rho)$ errors with $\omega$','Interpreter','latex')
legend('Interpreter','latex','Location',[0.4 0.52 0.1 0.25])



% Third graph

nexttile;hold on; grid on;
plot(log10(omega_seq), log10(1./m_h_eps_seq), '-', 'Color', blue, 'DisplayName', '$1/m_\omega^{h_\varepsilon}$', 'LineWidth', 4);
plot(log10(omega_seq), log10(1./m_2_seq), '-', 'Color',  purple, 'DisplayName', '$1/m_2$', 'LineWidth', 4);
plot(log10(omega_seq), log10(M_gamma_seq), '-.', 'Color',  'black', 'DisplayName', '$M_\omega^\gamma$', 'LineWidth', 4);
plot(log10(omega_seq), log10(sqrt(h_error_seq_square)), '-', 'Color', red, 'DisplayName', '$\|h - h_\varepsilon\|_{L^2(\rho)}$', 'LineWidth', 4);
plot(log10(omega_seq), log10(sqrt(g_error_seq_square)), '-.', 'Color', red, 'DisplayName', '$\|g - g_\varepsilon\|_{H^1_{\alpha}(\rho)}$', 'LineWidth', 4);

xlim([min(log10(omega_seq)), max(log10(omega_seq))])
legend('Interpreter','latex')
xlabel('$log_{10}(\omega)$','Interpreter','latex')
ylabel('$log_{10}$ of the corresponding values','Interpreter','latex')



set(gcf, 'position', [800, 800, 1200, 400])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',18)

print(fig, [fig_dir, 'conv_omega.pdf'], '-dpdf', '-r0')

%%
% %%
% figure; hold on; grid on;
% error_upper_bound_seq = 1./m_h_eps_seq.*(g_error_seq_square + M_gamma_seq.*h_error_seq_square);
% % loglog(omega_seq, upper_bound_errror, 'LineWidth', 4);
% plot(log10(omega_seq), log10(error_upper_bound_seq), '-', 'Color', 'green','DisplayName', 'error_upper_bound', 'LineWidth', 4);
% plot(log10(omega_seq), log10(1./m_h_eps_seq), '.', 'Color', 'green', 'DisplayName', 'coercivity_constant', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_regression_error_seq), '-', 'Color', 'black', 'DisplayName', 'regression error L^2(rho)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_prony_error_seq), '-',  'Color', 'blue', 'DisplayName', 'Prony error L^2(rho)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_regression_unif_error_seq),'-.', 'Color', 'black', 'DisplayName', 'regression error L^2(Lebesgue)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_prony_unif_error_seq), '-.',  'Color', 'blue', 'DisplayName', 'Prony error L^2(Lebesgue)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_regression_ill_error_seq),  'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(rho)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_regression_ill_unif_error_seq), '-.', 'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(Lebesgue)', 'LineWidth', 4);
% yline(log10(a))
% xlabel('log 10 omega')
% ylabel('log 10 L2 error')
% legend()
%
%
% %%
%
% figure; hold on; grid on;
% error_upper_bound_seq = 1./m_h_eps_seq.*(g_error_seq_square + M_gamma_seq.*h_error_seq_square);
% % loglog(omega_seq, upper_bound_errror, 'LineWidth', 4);
% plot(log10(omega_seq), log10(error_upper_bound_seq./theta_norm_omega), '-', 'Color', 'green','DisplayName', 'error_upper_bound', 'LineWidth', 4);
% plot(log10(omega_seq), log10(1./m_h_eps_seq./theta_norm_omega), '.', 'Color', 'green', 'DisplayName', 'coercivity_constant', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_regression_error_seq./theta_norm_omega), '-', 'Color', 'black', 'DisplayName', 'regression error L^2(rho)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_prony_error_seq./theta_norm_omega), '-',  'Color', 'blue', 'DisplayName', 'Prony error L^2(rho)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_regression_unif_error_seq./theta_norm_unif),'-.', 'Color', 'black', 'DisplayName', 'regression error L^2(Lebesgue)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_prony_unif_error_seq./theta_norm_unif), '-.',  'Color', 'blue', 'DisplayName', 'Prony error L^2(Lebesgue)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_regression_ill_error_seq./theta_norm_omega),  'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(rho)', 'LineWidth', 4);
% plot(log10(omega_seq), log10(gamma_regression_ill_unif_error_seq./theta_norm_unif), '-.', 'Color', 'red', 'DisplayName', 'ill-posed regression error L^2(Lebesgue)', 'LineWidth', 4);
% yline(log10(a))
% xlabel('log 10 omega')
% ylabel('log 10 L2 error')
% legend()
%
%
% %%
% figure;hold on;
% plot(log10(omega_seq), log10(h_error_seq_square./h_norm_omega))
% plot(log10(omega_seq), log10(g_error_seq_square./dg_norm_omega))
% legend('h rel err', 'dg rel err')
%
% %%
% figure;hold on;
% plot(log10(omega_seq),log10(loss_est_regression_seq))
% plot(log10(omega_seq),log10(loss_est_prony_seq))
% plot(log10(omega_seq),log10((gamma_regression_error_seq.*m_h_eps_seq).^2))
% plot(log10(omega_seq),log10((gamma_regression_error_seq.*M_h_eps_seq).^2))
% plot(log10(omega_seq),log10((gamma_prony_error_seq.*m_h_eps_seq).^2))
% plot(log10(omega_seq),log10((gamma_prony_error_seq.*M_h_eps_seq).^2))
% legend()
