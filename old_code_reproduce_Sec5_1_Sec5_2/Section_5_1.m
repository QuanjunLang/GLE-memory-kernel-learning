% Figure in Section 5.1:
% Trajectory and Noise
% Estimation of g and h
% Estimation of gamma, with F and without F
% Prediction error.


%% 
clc
close all
clear all


temp = matlab.desktop.editor.getActiveFilename;
parent_dir_temp = strsplit(temp, '/');
parent_dir = strjoin(parent_dir_temp(1:end-3), filesep);
fig_dir = [parent_dir, '/note/figure/'];


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
% F_type = 'steep_double_well';
% F_type = 'single_well';
F_type = '0';
switch F_type
    case '0'
        F = @(x) 0*x;
    case 'double_well'
        E0 = 0.02;
        F = @(x) -E0*x.*(x.^2 - 1);n 
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
traj_para.beta = 1;

M = 2^16;

figure;fplot(F)

tic;[v, dt, tgrid, T0, correlated_noise] = generate_trajectories(R, S, F, M, traj_para, random_seed);toc;
T = M*dt;
%% sparse and noisy observation of v
% obs_gap = 40;
% obs_len = floor(M/obs_gap/2);
% obs_std = 1e-4;
% VACF_obs_len = 30
obs_gap = 70;
obs_len = floor((M/obs_gap)/2);
obs_std = 1e-1;


tgrid_obs = tgrid(1:obs_gap:obs_len*obs_gap);
dt_obs = obs_gap *dt;
T_obs = obs_len*dt_obs;
v_obs = v(1:obs_gap:obs_len*obs_gap) + randn(1, obs_len)*obs_std;
M_obs = obs_len;


VACF_obs_len = 24;
%% Measure rho

omega = 0.05;
rho = @(x) exp(-2*omega*x); 


%% Figure of true kernel and trajectory
[fig, fig2] = plot_noise_traj_kernel(R, S, M, traj_para, correlated_noise, 28, v, v_obs, tgrid_obs, T0, fig_dir);
print(fig2,  [fig_dir, 'traj_noise.pdf'], '-dpdf', '-r0')
print(fig,  [fig_dir, 'kernel_spectrum.pdf'], '-dpdf', '-r0')

%% Generate VACF

corr_func = get_all_correlations(v, dt, F);
% RR = corr_func.h(1:L);
% g_dif_RR = corr_func.f(1:L);
% g_RR = corr_func.g(1:L);

corr_func_obs = get_all_correlations(v_obs, dt_obs, F);
% RR_obs = corr_func_obs.h(1:L_obs);
% g_dif_RR_obs = corr_func_obs.f(1:L_obs);
% g_RR_obs = corr_func_obs.g(1:L_obs);

%%



%% Prony method for VACF h
% Assume a finite observation of h,
% Use prony method to derive the analytical h


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


% a1 = sum(result_h.h(tgrid))*dt;
a2 = real(-sum(result_h.w./result_h.lam));
a1 = result_h.h(0);
comb_coef = a1^2/(a1^2 + a2^2);
% comb_coef is alpha_1 in the paper

%
Prony_g = result_g.g;
Prony_dg = result_g.dg;
Prony_h = result_h.h;



% error
h_error = sqrt(sum((corr_func.h - Prony_h(tgrid')).^2.*rho(tgrid'))*dt);
g_error = sqrt(sum((corr_func.g - Prony_g(tgrid')).^2.*rho(tgrid'))*dt);
dg_error = sqrt(sum((corr_func.dg - Prony_dg(tgrid')).^2.*rho(tgrid'))*dt);




%
h_error_paper = sum((corr_func.h - Prony_h(tgrid')).^2.*rho(tgrid'))*dt;
g_error_paper = comb_coef*sum((corr_func.g - Prony_g(tgrid')).^2.*rho(tgrid'))*dt + (1-comb_coef)*sum((corr_func.dg - Prony_dg(tgrid')).^2.*rho(tgrid'))*dt;



h_error_paper
g_error_paper

%% plot for Prony
% h
% figure; subplot(221);hold on;
% plot(tgrid_short, RR, 'LineWidth', 3)
% plot(tgrid_obs(1:VACF_obs_len), RR_obs(1:VACF_obs_len), 'o', 'LineWidth', 3)
% plot(tgrid_short, Prony_h(tgrid_short), 'LineWidth', 3)
% xlim([0, 20])
% legend('True VACF', 'Obs VACF', 'Est VACF')
% 
% 
% %g
% 
% subplot(222);hold on;
% plot(tgrid_short, g_RR, 'LineWidth', 4);
% plot(tgrid_short, Prony_g(tgrid_short), '.', 'Linewidth', 4);
% plot(tgrid_obs(1:VACF_gen_len), g_RR_obs, 'LineWidth', 4);
% % plot(tgrid_obs(1:VACF_gen_len), result., 'LineWidth', 4);
% xlim([0, 10])
% legend('g from RR (True)', 'g from Prony (The one we should use)', 'g derived from d RR obs', 'g estimted from obs traj')
% 
% %dg
% 
% subplot(223);hold on;
% plot(tgrid_short, dg_true, 'LineWidth', 4 );
% plot(tgrid_short, Prony_dg(tgrid_short), '.', 'LineWidth', 4 )
% plot(tgrid_obs(1:VACF_gen_len), gradient(g_RR_obs, dt_obs), 'LineWidth', 4 );
% 
% xlim([0, 20])
% legend('dg from RR (True)', 'dg from Prony h (The one we should use)', 'dg from RR obs')
% 
% 
% % modes
% subplot(224);hold on;
% scatter(real(result_h.r), imag(result_h.r), 130, '.'); hold on;circle(0, 0, 1)
% xrange = 1.4;
% xlim([-xrange, xrange]);ylim([-xrange, xrange])
% title('Prony modes')
% 
% 
% if F_type ~= '0'
%     figure;
%     
%     subplot(221);hold on;
%     plot(tgrid_short, g_dif_RR, 'LineWidth', 4 );
%     plot(tgrid_short, result_f.h(tgrid_short), 'o', 'LineWidth', 2 )
%     plot(tgrid_obs(1:VACF_gen_len), g_dif_RR_obs, 'LineWidth', 4 );
%     
%     xlim([0, 20])
%     legend('dg from RR (True)', 'dg from Prony h (The one we should use)', 'dg from RR obs')
%     
%     subplot(222);hold on;
%     plot(tgrid_short, gradient(g_dif_RR, dt), 'LineWidth', 4 );
%     plot(tgrid_short, result_f.g(tgrid_short), 'o', 'LineWidth', 2 )
%     % plot(tgrid_obs(1:VACF_gen_len), g_dif_RR_obs, 'LineWidth', 4 );
%     
%     xlim([0, 20])
%     legend('dg from RR (True)', 'dg from Prony h (The one we should use)', 'dg from RR obs')
%     
%     subplot(223);hold on;
%     plot(tgrid_short, gradient(gradient(RR, dt), dt), 'LineWidth', 4 );
%     plot(tgrid_short, result_h.dg(tgrid_short), 'o', 'LineWidth', 2 )
%     % plot(tgrid_obs(1:VACF_gen_len), g_dif_RR_obs, 'LineWidth', 4 );
%     
%     xlim([0, 20])
%     legend('dg from RR (True)', 'dg from Prony h (The one we should use)', 'dg from RR obs')
%     
% end
% 

%%
fig3 = figure;
t = tiledlayout(1, 3, "TileSpacing","compact");
t.TileSpacing = 'compact';
t.Padding = 'none';



nexttile;hold on;grid on;
plot(tgrid, corr_func.h, 'LineWidth', 3)
plot(tgrid, Prony_h(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.h(1:VACF_obs_len), '-o', 'LineWidth', 2)
xlim([0, 20])
title("h: Auto covariance of v");
legend('True h', 'Estimated h', 'Observed h')
xlabel('time t')

nexttile;hold on;grid on;
plot(tgrid, corr_func.g, 'LineWidth', 3)
plot(tgrid, Prony_g(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.g(1:VACF_obs_len), '-o', 'LineWidth', 2)
xlim([0, 20])
title("g: g = -h'");
legend('True g', 'Estimated g', 'Observed g')
xlabel('time t')


nexttile;hold on;grid on;
plot(tgrid, corr_func.dg, 'LineWidth', 3);
plot(tgrid, Prony_dg(tgrid), '-.', 'LineWidth', 4)
plot(tgrid_obs(1:VACF_obs_len), corr_func_obs.dg(1:VACF_obs_len), '-o', 'LineWidth', 2)

xlim([0, 20])
title("g': g' = -h''");
legend("True g'", "Estimated g'", "Observed g'", 'location', 'northeast')
xlabel('time t')





set(gcf, 'position', [300, 300, 800, 220])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)


print(fig3, [fig_dir, 'h_g_est.pdf'], '-dpdf', '-r0')


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
% regression_para.lb = 0;
% regression_para.rb = 30;
% regression_para.knot_num = 30;
% regression_para.deg = 3;
% regression_para.free_bdry = 1;
% regression_para.N = 2001;
% regression_para.omega = omega;
% regression_para.reg_method = 'LS';
% regression_para.h_conv_phi_dx_method = 'direct';
% regression_para.basis_type = 'spline';
% regression_para.lam = result_gamma.lam_seq;
% 
% h_est = @(x) real(Prony_h(x));
% dg_est = @(x) real(Prony_dg(x));
% g_est = @(x) real(Prony_g(x));
% 
% [theta, result_details_regression] = get_theta_regression(regression_para, h_est, dg_est);
% theta_ill_posed = get_theta_regression_ill_posed(regression_para, h_est, g_est);
% 
% xgrid_regression = linspace(regression_para.lb, regression_para.rb, regression_para.N)';
% dx_regression = xgrid_regression(2) - xgrid_regression(1);

%% Combined loss regression
 
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



% a1 = sum(result_h.h(tgrid))*dt;
% a2 = result_h.h(0);
% comb_coef = a1^2/(a1^2 + a2^2);
% comb_coef = 0.5;

regression_para.comb_coef = comb_coef;

h_est = @(x) real(Prony_h(x));
dg_est = @(x) real(Prony_dg(x));
g_est = @(x) real(Prony_g(x));

theta_derivative = get_theta_regression(regression_para, h_est, dg_est);
theta_ill_posed = get_theta_regression_ill_posed(regression_para, h_est, g_est);
[theta, result_details_regression] = get_theta_regression_combined(regression_para, h_est, dg_est, g_est);


xgrid_regression = linspace(regression_para.lb, regression_para.rb, regression_para.N)';
dx_regression = xgrid_regression(2) - xgrid_regression(1);

% figure;hold on;
% plot(xgrid_regression, theta(xgrid_regression), 'LineWidth', 4);
% plot(xgrid_regression, theta_ill_posed(xgrid_regression), 'LineWidth', 4);
% plot(xgrid_regression, theta_combined(xgrid_regression), 'LineWidth', 4);
% plot(xgrid_regression, R(xgrid_regression), 'LineWidth', 4);
% legend('Derivative', 'Spacial', 'Combined' ,'True')

% gamma_combined_error = sqrt(sum((theta_combined(tgrid_short') - R(tgrid_short')).^2.*rho(tgrid_short'))*dt);
% gamma_combined_error_2 = gamma_combined_error^2









%% Loss of the regression result, using estimated h and estimated g
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

% upbd = 1/m_h_eps*dg_error + M_gamma/m_h_eps*h_error;
% upbd


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

%% Strong traj compare
% figure; hold on;
% plot(tgrid_obs, pred_v_regression(1:obs_gap:obs_len*obs_gap));
% plot(tgrid_obs, pred_v_Prony(1:obs_gap:obs_len*obs_gap));
% plot(tgrid_obs, v(1:obs_gap:obs_len*obs_gap))
% plot(tgrid_obs, v_obs)
% legend('pred traj regression', 'pred traj prony', 'True traj', 'obs traj')
% 
% 
% traj_error_regression = sqrt(sum((pred_v_regression(1:obs_gap:obs_len*obs_gap) - v_obs).^2*dt_obs)/T_obs)
% traj_error_prony = sqrt(sum((pred_v_Prony(1:obs_gap:obs_len*obs_gap) - v_obs).^2*dt_obs)/T_obs)
% traj_error_true = sqrt(sum((v(1:obs_gap:obs_len*obs_gap) - v_obs).^2*dt_obs)/T_obs)
% 
% sqrt(sum((pred_v_Prony - v(1:M/2)).^2*dt)/T*2)
% sqrt(sum((pred_v_regression - v(1:M/2)).^2*dt)/T*2)
% 
% sqrt(sum((v(1:M/2)).^2*dt)/T*2)


%% Compare the vacf of trajectories generated by the estimated kernel
% basis_fourier_transform = spline_fourier_transform(regression_para);
% theta_Fourier_transform = @(y) 0*y;
% for j = 1:length(c_regression)
% theta_Fourier_transform = @(y) theta_Fourier_transform(y) + c_regression(j)*basis_fourier_transform{j}(y);
% end
% 
% %%
% theta_Fourier_transform_true = @(y) theta_Fourier_transform(y) + conj(theta_Fourier_transform(y));
% 





%% Coercivity analysis
save_coercivity = true;
new_coercivity;

%% error
gamma_regression_ill_posed_error = sqrt(sum((theta_ill_posed(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)
gamma_regression_derivative_error = sqrt(sum((theta_derivative(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)
gamma_prony_error = sqrt(sum((theta_Prony(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)
gamma_regression_combined_error = sqrt(sum((theta(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)



%% Plot estimation result
fig4 = figure;
hold on;
left_color = [.5 .5 0];
right_color = [0 .5 .5];
set(fig4,'defaultAxesColorOrder',[left_color; right_color]);
grey = 0.7*ones(1, 3);

yyaxis right
fill_rho_xgrid = [xgrid_regression; fliplr(xgrid_regression')'];
fill_rho_val = [rho(xgrid_regression); zeros(size(xgrid_regression))];
fill(fill_rho_xgrid, fill_rho_val, 'k', 'FaceAlpha', 0.15)
ylim([0, 2])
ylabel('measure \rho')

yyaxis left
grid on;
plot(xgrid_regression, R(xgrid_regression), 'LineWidth', 5, 'color', [0 0.4470 0.7410]);
plot(xgrid_regression, theta(xgrid_regression), '-.', 'LineWidth', 5, 'color', [0.8500 0.3250 0.0980]);
plot(xgrid_regression, theta_ill_posed(xgrid_regression), '-', 'LineWidth', 2, 'color', [0.4660    0.6740    0.1880]);
plot(xgrid_regression, theta_derivative(xgrid_regression), '-', 'LineWidth', 2, 'color',   [0.4940    0.1840    0.5560]);
plot(xgrid_regression(1:20:end), real(result_gamma.gamma(xgrid_regression(1:20:end))), '.', 'MarkerSize', 10, 'linewidth', 4, 'color', [0.9290 0.6940 0.1250]);
ylabel('memory kernel')
% h = get(gca,'Children');

% plot(xgrid_regression, theta_ill_posed(xgrid_regression), 'LineWidth', 4);

ylim([min(R(xgrid_regression))*1.5, max(R(xgrid_regression))*1.5])
legend()
legend('True kernel $\gamma$', ...
    ['$\theta, \|\theta -\gamma\|_{L^2(\rho)} = $', num2str(gamma_regression_combined_error)], ...
    ['$\theta_1, \|\theta_1 -\gamma\|_{L^2(\rho)} = $', num2str(gamma_regression_ill_posed_error)], ...
    ['$\theta_2, \|\theta_2 -\gamma\|_{L^2(\rho)} = $', num2str(gamma_regression_derivative_error)], ...
    ['$\theta_L, \|\theta_L -\gamma\|_{L^2(\rho)} = $', num2str(gamma_prony_error)], ...
    'measure $\rho$','Interpreter','latex')

title('Estimation result of the memory kernel')


ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = grey;

set(gcf, 'position', [300, 300, 700, 300])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)


print(fig4, [fig_dir, 'kernel_est.pdf'], '-dpdf', '-r0')

%%


reg_traj_para = traj_para;
reg_traj_para.method = 'reg';
reg_traj_para.loadON = 1;
theta_regression_fourier = result_details_regression.theta_Fourier_transform;
[predict_err_reg.traj, predict_err_reg.noise, predict_err_reg.h] = generate_traj_est_kernel(R, theta_regression_fourier, F, M, reg_traj_para, random_seed);


reg_traj_para = traj_para;
reg_traj_para.method = 'prony';
reg_traj_para.loadON = 1;
theta_prony_fourier = result_gamma.gamma_fourier;
[predict_err_prony.traj, predict_err_prony.noise, predict_err_prony.h] = generate_traj_est_kernel(R, theta_prony_fourier, F, M, reg_traj_para, random_seed);



%%

fig5 = figure;
hold on; grid on
% t = tiledlayout(1, 3, "TileSpacing","compact");
% t.TileSpacing = 'compact';
% t.Padding = 'none';



% nexttile;hold on;grid on;
plot(tgrid, corr_func.h, 'LineWidth', 3, 'LineWidth', 3, 'color', [0 0.4470 0.7410]);
plot(tgrid, predict_err_reg.h, '-.', 'LineWidth', 4, 'color', [0.8500 0.3250 0.0980]);
plot(tgrid(1:15:end), predict_err_prony.h(1:15:end), '.', 'MarkerSize', 8, 'linewidth', 4, 'color', [0.9290 0.6940 0.1250]);

xlim([0, 20])
title("Prediction error: VACF");
legend('True h', 'h: kernel estimated from regression', 'h: kernel estimated from Prony')
% xlabel('time t')


set(gcf, 'position', [300, 300, 300, 200])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)


print(fig5, [fig_dir, 'predict_error_h.pdf'], '-dpdf', '-r0')


%% Fourier transform of estimated memory kernels

fourier_grid = 0.01:0.01:8;

fig6 = figure;
hold on; grid on


plot(fourier_grid, S(fourier_grid), 'linewidth', 3, 'LineWidth', 3, 'LineWidth', 3, 'color', [0 0.4470 0.7410]);
plot(fourier_grid, result_details_regression.theta_Fourier_transform(fourier_grid), '-.', 'LineWidth', 4, 'color', [0.8500 0.3250 0.0980]);
plot(fourier_grid, result_gamma.gamma_fourier(fourier_grid), '.', 'MarkerSize', 8, 'linewidth', 4, 'color', [0.9290 0.6940 0.1250]);
legend('S', 'approx')

xlim([0, fourier_grid(end)])
title("Prediction error: spectral functions");
legend('True', 'Regression', 'Inv Laplace Transform')
% xlabel('time t')


set(gcf, 'position', [300, 300, 300, 200])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)


print(fig6, [fig_dir, 'predict_error_S.pdf'], '-dpdf', '-r0')

%%
figure;hold on;
plot(v);plot(predict_err_prony.traj);plot(predict_err_reg.traj)




