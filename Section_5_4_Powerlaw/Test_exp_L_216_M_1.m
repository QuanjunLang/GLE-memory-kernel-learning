
function result = Test_exp_L_216_M_1()

random_seed = 5;
rng(random_seed);
addPathAll;



%% Generate true memory kernel along with it s Fourier transform and Laplace transform
memoryInfo.kernel_type = 'Test';
% memoryInfo.kernel_type = 'RandomProny_p_2_q_2'; 
% 
% memoryInfo.kernel_type  = 'PowerLaw';
% memoryInfo.kernel_type  = 'Polynomial';
memoryInfo.plotON       = 1;

memoryInfo = generate_memory_kernel(memoryInfo);

%% Generate system info
sysInfo.m      = 1;
sysInfo.beta   = 1;
sysInfo.mu     = 0;

sysInfo.F_type = 'linear';      % linear, double_well
sysInfo.G_type = '0';

sysInfo.fig_dir     = fig_dir;
sysInfo.data_dir    = data_dir;

sysInfo = update_system_info(sysInfo, memoryInfo);  % Generate the true auto correlation when F is linear

%% Generate correlated noise and trajectories
trajInfo.loadON    = 1;
trajInfo.plotON    = 0;
trajInfo.wu        = 80*pi;       % Spectral cutoff
trajInfo.N         = 8000;        % Number of spectral points
trajInfo.L         = 2^16;
trajInfo.dt        = pi/trajInfo.wu;
trajInfo.M         = 1;

trajInfo = generate_trajectories(memoryInfo, sysInfo, trajInfo, random_seed);

sysInfo.dt      = trajInfo.dt;
sysInfo.tgrid   = trajInfo.tgrid;
%% sparse and noisy observation of v
obsInfo.plotON  = 1;
obsInfo.gap     = 10;
obsInfo.std     = 0.01;
obsInfo.len     = floor(trajInfo.L/obsInfo.gap);

obsInfo = generate_observed_trajectory(obsInfo, trajInfo);
%% Estimate correlation function from discrete trajectories
sysInfo.omega  = 0.1;                     % Define Measure rho
sysInfo.rho    = @(x) exp(-2*sysInfo.omega*x);


temporal_spacial_method = 'M_only';

[corr_func, corr_func_obs] = get_correlations_true_obs(sysInfo, trajInfo, obsInfo, temporal_spacial_method);



%% Interpolate using Prony method
% Assume a finite observation of h,
% Use prony method to derive the analytical h


VACF_obs_len = 50;
% options for the velocity autocorrelation
PronyInfo_h.prony_p             = 10;
PronyInfo_h.prony_N             = VACF_obs_len;
PronyInfo_h.obs_dx              = obsInfo.dt_obs;
PronyInfo_h.obs_h_grid          = corr_func_obs.h(1:VACF_obs_len);
PronyInfo_h.rho                 = sysInfo.rho;
PronyInfo_h.tgrid_obs           = obsInfo.tgrid_obs;
PronyInfo_h.polycoef_method     = 'MP';
PronyInfo_h.weight_method       = 'LS_h0_new';
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

%% Estimate Gamma using Laplace transform

switch sysInfo.mu
    case 0
        result_gamma = get_gamma_from_prony_h(PronyInfo_h, PronyDetails.result_h);
        theta_Prony = @(x) real(result_gamma.gamma(x));
    otherwise
        result_gamma = get_gamma_from_prony_h_esternal_force(PronyInfo_h, PronyInfo_g, PronyDetails.result_h, PronyDetails.result_phi);
        theta_Prony = @(x) real(result_gamma.gamma(x));
end

figure;
subplot(121);hold on;

R = memoryInfo.kernel;
gamma = R;
xgrid = obsInfo.tgrid_obs;
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


%% Combined loss regression

regression_para.lb          = 0;
regression_para.rb          = 30;
regression_para.knot_num    = 50;
regression_para.deg         = 3;
regression_para.free_bdry   = 1;
regression_para.N           = 2001;
regression_para.omega       = sysInfo.omega;
regression_para.reg_method  = 'LS';

regression_para.convolution_method  = 'direct';
% regression_para.basis_type          = 'prony';
regression_para.basis_type          = 'spline';
regression_para.lam                 = result_gamma.lam_seq;



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


%% Loss of the regression result, using estimated h and estimated g
regression_debug


%% Coercivity analysis
save_coercivity = false;

omega = sysInfo.omega;
R_lap = memoryInfo.kernel_Laplace;
% R_lap = @(s) temp_func(s);
comb_coef = EstResult.comb_coef;
result_h = PronyDetails.result_h;
h_error_paper = EstResult.error.h^2;
g_error_paper = EstResult.error.g^2;

new_coercivity;


%% Plot estimation result
tgrid   = sysInfo.tgrid;
dt      = sysInfo.dt;

gamma_regression_1_error    = sqrt(sum((theta_1(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
gamma_regression_2_error    = sqrt(sum((theta_2(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
gamma_regression_0_error    = sqrt(sum((theta_0(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
gamma_prony_error           = sqrt(sum((theta_Prony(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);


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
plot(xgrid_regression(1:20:end), real(result_gamma.gamma(xgrid_regression(1:20:end))), '.', 'MarkerSize', 10, 'linewidth', 4, 'color', [0.9290 0.6940 0.1250]);
ylabel('memory kernel')
% h = get(gca,'Children');

% plot(xgrid_regression, theta_ill_posed(xgrid_regression), 'LineWidth', 4);

ylim([min(R(xgrid_regression))*1.5, max(R(xgrid_regression))*1.5])
legend()
legend('True kernel $\gamma$', ...
    ['$\theta, \|\theta -\gamma\|_{L^2(\rho)} = $', num2str(gamma_regression_0_error)], ...
    ['$\theta_1, \|\theta_1 -\gamma\|_{L^2(\rho)} = $', num2str(gamma_regression_1_error)], ...
    ['$\theta_2, \|\theta_2 -\gamma\|_{L^2(\rho)} = $', num2str(gamma_regression_2_error)], ...
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


% print(fig4, [fig_dir, 'kernel_est.pdf'], '-dpdf', '-r0')


%%
g_est;
h_est;

% Inputs: function handles h_est(t), g_est(t)

% Time discretization
T = 200;                 % Final time
N = 2^12;               % Number of time points (power of 2 for FFT)
dt = T / N;
tgrid = (0:N-1)' * dt;  % Column vector

% Evaluate functions
h_vals = h_est(tgrid);
g_vals = g_est(tgrid);

% Zero-pad to length 2N to minimize circular convolution artifacts
L = 2 * N;
h_pad = [h_vals; zeros(N,1)];
g_pad = [g_vals; zeros(N,1)];

% FFT
H_fft = fft(h_pad);
G_fft = fft(g_pad);

% Tikhonov regularization in Fourier domain
lambda = 1e-2;  % Regularization parameter
Gamma_fft = conj(H_fft) .* G_fft ./ (abs(H_fft).^2 + lambda);

% Inverse FFT to time domain and correct scaling
gamma_rec = real(ifft(Gamma_fft)) / dt;

% Truncate to original interval
gamma_vals = gamma_rec(1:N);

% Optionally return as function handle
theta_F = @(t_query) interp1(tgrid, gamma_vals, t_query, 'linear', 0);

% Plot result
figure;hold on;
plot(tgrid, gamma_vals, 'LineWidth', 2);
plot(tgrid, theta_0(tgrid), 'LineWidth', 2);
xlabel('t'); ylabel('\gamma(t)');
title('Recovered \gamma(t) using Tikhonov–Fourier Deconvolution');
grid on;
xlim([0, 10])


error_F = sqrt(sum((theta_F(tgrid) - memoryInfo.kernel(tgrid)).^2 .*rho(tgrid)));



%% 


close all

result.h_est = h_est;
result.g_est = g_est;
result.dg_est = dg_est;

tgrid = sysInfo.tgrid;
dt = sysInfo.dt;

result.theta_0 = theta_0;
result.theta_1 = theta_1;
result.theta_2 = theta_2;
result.theta_L = theta_Prony;
result.theta_F = theta_F;
result.theta_T = memoryInfo.kernel;

result.rho      = sysInfo.rho;
result.tgrid    = sysInfo.tgrid;

result.error_0 = gamma_regression_0_error;
result.error_1 = gamma_regression_1_error;
result.error_2 = gamma_regression_2_error;
result.error_L = gamma_prony_error;
result.error_F = error_F;

result.h_true = sysInfo.sym.h_true;
result.g_true = sysInfo.sym.g_true;
result.dg_true = sysInfo.sym.dg_true;

result.h_error = sqrt(sum((result.h_est(tgrid) - result.h_true(tgrid)).^2 .*rho(tgrid) *dt));
result.g_error = sqrt((1-comb_coef)*sum((result.g_est(tgrid) - result.g_true(tgrid)).^2.*rho(tgrid))*dt + comb_coef*sum((result.dg_est(tgrid) - result.dg_true(tgrid)).^2.*rho(tgrid))*dt);

result.comb_coef = comb_coef;
result.upbd = upbd_iprv;



end