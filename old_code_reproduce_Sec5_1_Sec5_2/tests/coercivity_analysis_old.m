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




a1 = sum(result_h.h(tgrid))*dt;
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