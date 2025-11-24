NN = 2000;
h_true_short = RR(1:NN);
dg_true_short = dg_true(1:NN);

tgrid_short = tgrid(1:NN)';
rho_true_short = rho(tgrid_short);

% Loss of the regression result, using estimated h and estimated g
% loss_est_regression  = sum((h_conv_phi_dx*c + dg_grid_regression).^2.*rho_grid_regression)*dx

omega_seq = linspace(0.001, 1, 100);
L2error_seq = zeros(100, 1);
L2error_coer_seq = zeros(100, 1);
loss_seq = zeros(100, 1);
for i = 1:100
    i
    omega = omega_seq(i);
    rho = @(x) exp(-omega*x);        % omega

% hgrid_regression = result_prony.h(xgrid_regression);
rho_grid_regression = rho(xgrid_regression);
    theta = regression_omega(omega, result_prony);
    
    
    conv_temp = zeros(N, 1);
    conv_temp(1) = 0;
    phi = theta(xgrid_regression);
    for n = 2:N
        val = (hgrid_regression(1)*phi(n) + hgrid_regression(n)*phi(1))/2;
        val = val + sum(hgrid_regression(2:n-1) .* phi(n-1:-1:2));
        conv_temp(n) = val*dx;
    end
    
    R_conv_h_method1 = conv_temp;
    R_conv_h_dx_method1 = gradient(R_conv_h_method1, dx);
    
    loss_true = sum((R_conv_h_dx_method1 + dg_grid_regression).^2.*rho_grid_regression)*dx;
    
    figure;hold on;plot(xgrid_regression,R_conv_h_dx_method1);plot(xgrid_regression, -dg_grid_regression, 'o')
    figure;hold on;plot(xgrid_regression, R_conv_h_method1);plot(xgrid_regression, -Prony_g(xgrid_regression))
    figure;hold on;plot(xgrid_regression, R_conv_h_dx_method1);plot(tgrid_short, -dg_true_short)
    % figure;
    
    
    % Secondly, the bound for laplace transform of h_eps
    tau_grid_real = linspace(-200, 200, 1000);
    tau_grid = omega + 1i*tau_grid_real;
    z_h_eps_z = abs(tau_grid.*result_prony.h_lap(tau_grid));
    
    m_h_eps = min(z_h_eps_z);
    M_h_eps = max(z_h_eps_z);
    
    gamma_regression_error = sqrt(sum((theta(tgrid') - R(tgrid')).^2.*rho(tgrid'))*dt);
    L2error_seq(i) = gamma_regression_error;
    L2error_coer_seq(i) = (m_h_eps * gamma_regression_error)^2;
    loss_seq(i) = loss_true;
end
%%
figure;hold on;
plot(omega_seq, loss_seq);
plot(omega_seq, L2error_coer_seq, 'o');

figure;
plot(omega_seq, L2error_seq, 'o');
%% Loss of the regression result, using estimated h and estimated g

conv_temp = zeros(N, 1);
conv_temp(1) = 0;
phi = theta(xgrid_regression);
ggrid_regression = Prony_g(xgrid_regression);

for n = 2:N
    val = (ggrid_regression(1)*phi(n) + ggrid_regression(n)*phi(1))/2;
    val = val + sum(ggrid_regression(2:n-1) .* phi(n-1:-1:2));
    conv_temp(n) = val*dx;
end
conv_temp = conv_temp + hgrid_regression(1)*phi;

% R_conv_h_method1 = conv_temp;
R_conv_h_dx_method2 = conv_temp;

loss_true = sum((R_conv_h_dx_method2 + dg_grid_regression).^2.*rho_grid_regression)*dx

figure;hold on;plot(xgrid_regression,R_conv_h_dx_method2);plot(xgrid_regression, -dg_grid_regression, 'o')
% figure;hold on;plot(xgrid_regression, R_conv_h_method1);plot(xgrid_regression, -Prony_g(xgrid_regression))
figure;hold on;plot(xgrid_regression, R_conv_h_dx_method2);plot(tgrid_short, -dg_true_short)
figure;hold on;plot(xgrid_regression, theta(xgrid_regression));plot(xgrid_regression, R(xgrid_regression))




%%
dif = @(x) R(x) - theta(x);

figure;plot(xgrid_regression, dif(xgrid_regression))
sqrt(sum(dif(xgrid_regression).^2.*rho(xgrid_regression))*dx)











%%

conv_temp = zeros(NN, 1);
conv_temp(1) = 0;
phi = R(tgrid_short);
for n = 2:NN
    val = (h_true_short(1)*phi(n) + h_true_short(n)*phi(1))/2;
    val = val + sum(h_true_short(2:n-1) .* phi(n-1:-1:2));
    conv_temp(n) = val*dt;
end

R_conv_h = conv_temp;
R_conv_h_dx = gradient(R_conv_h, dt);

loss_true = sum((R_conv_h_dx + dg_true_short).^2.*rho_true_short)*dt



% Loss of the regression result, using estimated h and estimated g
%
% conv_temp = zeros(N, 1);
% conv_temp(1) = 0;
% phi = theta(xgrid_regression);
% for n = 2:N
%     val = (hgrid_regression(1)*phi(n) + hgrid_regression(n)*phi(1))/2;
%     val = val + sum(hgrid_regression(2:n-1) .* phi(n-1:-1:2));
%     conv_temp(n) = val*dx;
% end
% theta_conv_h = conv_temp;
% theta_conv_h_dx = gradient(theta_conv_h, dx);

% loss_true_regression = sum((theta_conv_h_dx + dg_true).^2.*rho(tgrid'))*dt



% Loss of the regression result, using estimated h and estimated g
% loss_est_prony  = sum((theta_Prony(xgrid_regression) + dg_grid_regression).^2.*rho_grid_regression)*dx
% loss_true_prony = sum((theta_Prony(tgrid') + dg_true).^2.*rho(tgrid'))*dt