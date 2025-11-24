function result_all = ensemble_algorithm(sys_para, traj_para, obs_para, Prony_I, Prony_I_f, regression_para)



%% Generate correlated noise and trajectories
% F_type = 'steep_double_well';
% F_type = 'single_well';
F_type = 'very_steep_double_well';

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
    case 'very_steep_double_well'
        E0 = 1;
        F = @(x) -E0*x.*(x.^2 - 4);
    case 'very_steep_and_far_double_well'
        E0 = 1;
        F = @(x) -E0*x.*(x.^2 - 16);
end




G_type = '0';
switch G_type
    case '0'
        G = @(x) 0*x;
    case 'Duffing'
        G_gamma = 0.1;
        G_omega = 1;
        G = @(x) G_gamma*cos(G_omega*x);
    case 'Strong_Duffing'
        G_gamma = 1;
        G_omega = 1;
        G = @(x) G_gamma*cos(G_omega*x);
    case 'Extra_Strong_Duffing'
        G_gamma = 10;
        G_omega = 10;
        G = @(x) G_gamma*cos(G_omega*x);
    case 'Extra_Strong_slow_Duffing'
        G_gamma = 10;
        G_omega = 1;
        G = @(x) G_gamma*cos(G_omega*x);
end



traj_para.traj_num  = 1;             % available noise traj length: 50, 2000, 5000
traj_para.F_type    = F_type;
traj_para.loadON    = 0;
traj_para.saveON    = 0;
traj_para.wu        = 80*pi;
traj_para.N         = 8000;
traj_para.beta      = 0.3;
traj_para.G         = G;
traj_para.G_type    = G_type;
traj_para.F         = F;

M = 2^12; % or 2^12, 2^14

if traj_para.traj_num == 1
    M = 2^12;
    % M = 2^16;
    % [all_v, dt, tgrid, T0, all_correlated_noise] = generate_trajectories(R, S, F, M, traj_para, random_seed);
    % T = M*dt;
    % all_v = all_v';

    corr_method = 'single_traj';

else


    corr_method = 'mixed';
end

% corr_method = 'monte_carlo';
[all_v, dt, tgrid, T0, all_correlated_noise] = generate_trajectories_ensemble(R, S, F, M, traj_para, random_seed);
T = M*dt;

%% Plot the trajectory
plot_2d_ON = false;
plot_1d_ON = false;

plot_motion_1d_2d

%% sparse and noisy observation of v
obs_gap = 30;
obs_len = floor((M/obs_gap));
obs_std = 1e-3;

tgrid_obs = tgrid(1:obs_gap:obs_len*obs_gap);
dt_obs = obs_gap *dt;
T_obs = obs_len*dt_obs;
all_v_obs = all_v(1:obs_gap:obs_len*obs_gap, :) + randn(obs_len, traj_para.traj_num)*obs_std;
M_obs = obs_len;


%% correlation functions


[corr_func, corr_func_obs] = get_all_correlations_all_method(all_v, all_v_obs, dt, dt_obs, F, corr_method, R, tgrid);

%% Prony method for VACF h
% Assume a finite observation of h,
% Use prony method to derive the analytical h

VACF_obs_len = 70;
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

Prony_I_f = Prony_I;
Prony_I_f.weight_method = 'LS';
Prony_I_f.obs_h_grid = corr_func_obs.f(1:VACF_obs_len);
Prony_I_f.prony_p = 15;
Prony_I_f.drop_0 = 0;

[result_h, result_g, result_f] = estimate_h_g_prony(Prony_I, Prony_I_f, F_type);


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
f_error = sqrt(sum((corr_func.f - result_f.h(tgrid')).^2.*rho(tgrid'))*dt);
g_error = sqrt(sum((corr_func.g - Prony_g(tgrid')).^2.*rho(tgrid'))*dt);
dg_error = sqrt(sum((corr_func.dg - Prony_dg(tgrid')).^2.*rho(tgrid'))*dt);


h_error_paper = sum((corr_func.h - Prony_h(tgrid')).^2.*rho(tgrid'))*dt
g_error_paper = comb_coef*sum((corr_func.g - Prony_g(tgrid')).^2.*rho(tgrid'))*dt + (1-comb_coef)*sum((corr_func.dg - Prony_dg(tgrid')).^2.*rho(tgrid'))*dt


h_error
f_error
g_error
dg_error


%%
new_fig_3_h_g_dg_est

%% Estimate Gamma using Laplace transform
switch F_type
    case '0'
        result_gamma = get_gamma_from_prony_h(Prony_I, result_h);
        theta_Prony = @(x) real(result_gamma.gamma(x));
    otherwise
        result_gamma = get_gamma_from_prony_h_esternal_force(Prony_I, Prony_I_f, result_h, result_f);
        theta_Prony = @(x) real(result_gamma.gamma(x));
end

%%
 
regression_para.lb = 0;
regression_para.rb = 50;
regression_para.knot_num = 30;
regression_para.deg = 3;
regression_para.free_bdry = 1;
regression_para.N = 2001;
regression_para.omega = omega;
regression_para.reg_method = 'RKHS';
regression_para.h_conv_phi_dx_method = 'direct';
regression_para.basis_type = 'spline';      % or 'prony'%
regression_para.lam = result_gamma.lam_seq;

% a1 = sum(result_h.h(tgrid))*dt;
% a2 = result_h.h(0);
% comb_coef = a1^2/(a1^2 + a2^2);

regression_para.comb_coef = comb_coef;



h_est = @(x) real(Prony_h(x));
dg_est = @(x) real(Prony_dg(x));
g_est = @(x) real(Prony_g(x));


theta_derivative = get_theta_regression(regression_para, h_est, dg_est);
theta_ill_posed = get_theta_regression_ill_posed(regression_para, h_est, g_est);
[theta, result_details_regression] = get_theta_regression_combined(regression_para, h_est, dg_est, g_est);


xgrid_regression = linspace(regression_para.lb, regression_para.rb, regression_para.N)';
dx_regression = xgrid_regression(2) - xgrid_regression(1);

%% error

gamma_regression_ill_posed_error = sqrt(sum((theta_ill_posed(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)
gamma_regression_derivative_error = sqrt(sum((theta_derivative(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)
gamma_prony_error = sqrt(sum((theta_Prony(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)
gamma_regression_combined_error = sqrt(sum((theta(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt)


%%
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


% print(fig4, [fig_dir, 'kernel_est.pdf'], '-dpdf', '-r0')

%%
save_coercivity = false;
new_coercivity;

end