basis_mat       = result_details_regression_0.basis_mat;
xgrid           = result_details_regression_0.xgrid;
h_conv_phi      = result_details_regression_0.h_conv_phi;
h_conv_phi_dx   = result_details_regression_0.h_conv_phi_dx;


basis = basis_mat{1};

gamma_true = memoryInfo.kernel;

gamma_true_grid = gamma_true(xgrid);



c_oracal = basis\ gamma_true_grid;
c = result_details_regression_0.c;
A = result_details_regression_0.A;
b = result_details_regression_0.b;




h_conv_gamma_est = h_conv_phi*c;
h_conv_gamma_oracal = h_conv_phi*c_oracal;

g_Prony = g_est(xgrid);

h_conv_gamma_dx_est = h_conv_phi_dx*c;
h_conv_gamma_dx_oracal = h_conv_phi_dx*c_oracal;

dg_Prony = dg_est(xgrid);


%
figure;

subplot(1,2,1); hold on;

plot(xgrid, h_conv_gamma_est,  '-.', 'LineWidth', 7, 'DisplayName', 'Estimated h * \gamma');
plot(xgrid, h_conv_gamma_oracal,  'LineWidth', 4, 'DisplayName', 'Oracle h * \gamma');
if isfield(sysInfo, 'sym')
    g_true = sysInfo.sym.g_true(xgrid);
    plot(xgrid, g_true,               'LineWidth', 4, 'DisplayName', 'True g');
end
plot(xgrid, g_Prony,              'LineWidth', 4, 'DisplayName', 'Prony g');

legend('Location', 'best');
xlabel('x');
ylabel('Value');
title('Convolution Comparison');
grid on;

subplot(1,2,2); hold on;

plot(xgrid, h_conv_gamma_dx_est,  '-.',   'LineWidth',7 , 'DisplayName', 'Estimated d h * \gamma');
plot(xgrid, h_conv_gamma_dx_oracal, 'LineWidth', 4, 'DisplayName', 'Oracle d h * \gamma');
if isfield(sysInfo, 'sym')
    dg_true = sysInfo.sym.dg_true(xgrid);
    plot(xgrid, dg_true,                'LineWidth', 4, 'DisplayName', 'True dg');
end
plot(xgrid, dg_Prony,               'LineWidth', 4, 'DisplayName', 'Prony dg');

legend('Location', 'best');
xlabel('x');
ylabel('Value');
title('Derivative Convolution Comparison');
grid on;



m = sysInfo.m;
rho = sysInfo.rho;

hgrid_regression = h_est(xgrid_regression);
theta_grid_regression = theta_0(xgrid_regression);
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



figure;

% Left subplot: h convolution
subplot(1,2,1); hold on; grid on;
plot(h_conv_R, '-', 'LineWidth', 2, 'DisplayName', 'True h * R');
plot(h_conv_theta, '--', 'LineWidth', 2, 'DisplayName', 'Estimated h * θ');
plot(g_grid_regression, ':', 'LineWidth', 2.5, 'DisplayName', 'g');

legend('Location', 'best');
xlabel('Time');
ylabel('Function Value');
title('Convolution of h');
set(gca, 'FontSize', 12);

% Right subplot: dh convolution
subplot(1,2,2); hold on; grid on;
plot(h_conv_R_dx, '-', 'LineWidth', 2, 'DisplayName', 'True h'' * R');
plot(h_conv_theta_dx, '--', 'LineWidth', 2, 'DisplayName', 'Estimated h'' * θ');
plot(dg_grid_regression, ':', 'LineWidth', 2.5, 'DisplayName', 'g''');

legend('Location', 'best');
xlabel('Time');
ylabel('Function Value');
title('Convolution of h''');
set(gca, 'FontSize', 12);


a = 1;

