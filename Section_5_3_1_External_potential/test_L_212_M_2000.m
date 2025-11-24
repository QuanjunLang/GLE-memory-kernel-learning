function result = test_L_212_M_2000()

%% 
% % clc
% close all
% clear all
random_seed = 20;
rng(random_seed);
addPathAll;



%% Generate true memory kernel along with it s Fourier transform and Laplace transform
% memoryInfo.kernel_type = 'Test';
% memoryInfo.kernel_type = 'RandomProny_p_1_q_1'; 

memoryInfo.kernel_type  = 'PowerLaw';
% memoryInfo.kernel_type  = 'Polynomial';
memoryInfo.plotON       = 0;

memoryInfo = generate_memory_kernel(memoryInfo);

%% Generate system info
sysInfo.m      = 0.3;
sysInfo.beta   = 0.5;
sysInfo.mu     = 1;

sysInfo.F_type = 'very_steep_double_well';      % linear, double_well
sysInfo.G_type = '0';

sysInfo.fig_dir     = fig_dir;
sysInfo.data_dir    = data_dir;

sysInfo = update_system_info(sysInfo, memoryInfo);  % Generate the true auto correlation when F is linear

%% Generate correlated noise and trajectories
trajInfo.loadON    = 1;
trajInfo.plotON    = 0;
trajInfo.wu        = 80*pi;       % Spectral cutoff
trajInfo.N         = 8000;        % Number of spectral points
trajInfo.L         = 2^12;
trajInfo.dt        = pi/trajInfo.wu;
trajInfo.M         = 2000;

trajInfo = generate_trajectories(memoryInfo, sysInfo, trajInfo, random_seed);

sysInfo.dt      = trajInfo.dt;
sysInfo.tgrid   = trajInfo.tgrid;
%% sparse and noisy observation of v
obsInfo.plotON  = 1;
obsInfo.gap     = 10;
obsInfo.std     = 1e-2;
obsInfo.len     = floor(trajInfo.L/obsInfo.gap);

obsInfo = generate_observed_trajectory(obsInfo, trajInfo);
%% Estimate correlation function from discrete trajectories
sysInfo.omega  = 0.05;                     % Define Measure rho
sysInfo.rho    = @(x) exp(-2*sysInfo.omega*x);
rho = sysInfo.rho;

temporal_spacial_method = 'M_only';

[corr_func, corr_func_obs] = get_correlations_true_obs(sysInfo, trajInfo, obsInfo, temporal_spacial_method);



%% Interpolate using Prony method
% Assume a finite observation of h,
% Use prony method to derive the analytical h


VACF_obs_len = 30;
% options for the velocity autocorrelation
PronyInfo_h.prony_p             = 10;
PronyInfo_h.prony_N             = VACF_obs_len;
PronyInfo_h.obs_dx              = obsInfo.dt_obs;
PronyInfo_h.obs_h_grid          = corr_func_obs.h(1:VACF_obs_len);
PronyInfo_h.rho                 = sysInfo.rho;
PronyInfo_h.tgrid_obs           = obsInfo.tgrid_obs;
PronyInfo_h.polycoef_method     = 'MP';
PronyInfo_h.weight_method       = 'LS';
PronyInfo_h.root_normalization  = 1;
PronyInfo_h.lambda_augmentation = 1;
PronyInfo_h.drop_0              = 0;

% options for the forcing autocorrelation

PronyInfo_g                     = PronyInfo_h;
PronyInfo_g.weight_method       = 'LS';
PronyInfo_g.obs_h_grid          = corr_func_obs.f(1:VACF_obs_len);
PronyInfo_g.prony_p             = PronyInfo_h.prony_p;
PronyInfo_g.drop_0              = 0;

% Estimate using prony method
[EstResult, PronyDetails] = estimate_h_g_prony(PronyInfo_h, PronyInfo_g, sysInfo, corr_func, corr_func_obs);

Prony_g     = EstResult.g;
Prony_dg    = EstResult.dg;
Prony_h     = EstResult.h;

%% Combined loss regression
R = memoryInfo.kernel;

regression_para.lb          = 0;
regression_para.rb          = 30;
regression_para.knot_num    = 50;
regression_para.deg         = 3;
regression_para.free_bdry   = 1;
regression_para.N           = 2001;
regression_para.omega       = sysInfo.omega;
regression_para.reg_method  = 'RKHS';

regression_para.convolution_method  = 'direct';
% regression_para.basis_type          = 'prony';
regression_para.basis_type          = 'spline';



h_est   = @(x) real(Prony_h(x));
dg_est  = @(x) real(Prony_dg(x));
g_est   = @(x) real(Prony_g(x));



% EstResult.comb_coef = 0.8;

[theta_1, result_details_regression_1]  = get_theta_regression_combined(regression_para, h_est, dg_est, g_est, 0);
[theta_2, result_details_regression_2]  = get_theta_regression_combined(regression_para, h_est, dg_est, g_est, 1);
[theta_0, result_details_regression_0]  = get_theta_regression_combined(regression_para, h_est, dg_est, g_est, EstResult.comb_coef);


xgrid_regression = linspace(regression_para.lb, regression_para.rb, regression_para.N)';
dx_regression = xgrid_regression(2) - xgrid_regression(1);

figure;hold on;
plot(xgrid_regression, theta_2(xgrid_regression), '--', 'LineWidth', 5);
plot(xgrid_regression, theta_1(xgrid_regression), 'LineWidth', 4);
plot(xgrid_regression, theta_0(xgrid_regression), 'LineWidth', 4);
plot(xgrid_regression, R(xgrid_regression), 'LineWidth', 4);
legend('Derivative', 'Spacial', 'Combined' ,'True')




%% Plot estimation result
tgrid   = sysInfo.tgrid;
dt      = sysInfo.dt;

gamma_regression_1_error    = sqrt(sum((theta_1(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
gamma_regression_2_error    = sqrt(sum((theta_2(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
gamma_regression_0_error    = sqrt(sum((theta_0(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
% gamma_prony_error           = sqrt(sum((theta_Prony(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);


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
plot(xgrid_regression, theta_0(xgrid_regression), '-.', 'LineWidth', 5, 'color', [0.8500 0.3250 0.0980]);
plot(xgrid_regression, theta_1(xgrid_regression), '-', 'LineWidth', 2, 'color', [0.4660    0.6740    0.1880]);
plot(xgrid_regression, theta_2(xgrid_regression), '-', 'LineWidth', 2, 'color',   [0.4940    0.1840    0.5560]);
% plot(xgrid_regression(1:20:end), real(result_gamma.gamma(xgrid_regression(1:20:end))), '.', 'MarkerSize', 10, 'linewidth', 4, 'color', [0.9290 0.6940 0.1250]);
ylabel('memory kernel')
% h = get(gca,'Children');

% plot(xgrid_regression, theta_ill_posed(xgrid_regression), 'LineWidth', 4);

ylim([min(R(xgrid_regression))*1.5, max(R(xgrid_regression))*1.5])
legend()
legend('True kernel $\gamma$', ...
    ['$\theta, \|\theta -\gamma\|_{L^2(\rho)} = $', num2str(gamma_regression_0_error)], ...
    ['$\theta_1, \|\theta_1 -\gamma\|_{L^2(\rho)} = $', num2str(gamma_regression_1_error)], ...
    ['$\theta_2, \|\theta_2 -\gamma\|_{L^2(\rho)} = $', num2str(gamma_regression_2_error)],'measure $\rho$','Interpreter','latex')
    % ['$\theta_L, \|\theta_L -\gamma\|_{L^2(\rho)} = $', num2str(gamma_prony_error)], ...
    

title('Estimation result of the memory kernel')


ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = grey;

set(gcf, 'position', [300, 300, 700, 300])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(findall(gcf,'-property','FontSize'),'FontSize',14)




%%

% result.w_seq            = memoryInfo.w_seq;
% result.lambda_seq       = memoryInfo.lambda_seq;
result.theta_true = memoryInfo.kernel;
result.omega = sysInfo.omega;
result.theta = theta_0;
result.theta_error = gamma_regression_0_error;
result.comb_coef = EstResult.comb_coef;


end