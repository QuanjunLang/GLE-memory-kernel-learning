function [m, M] = get_m_M(omega, h)


% Firstly, the bound for laplace transform of gamma
tau_grid_real = linspace(-20, 20, 100000);
% tau_grid_real = 10.^linspace(-20, 20, 100000);
tau_grid = omega + 1i*tau_grid_real;
z_gamma_z = abs(tau_grid.*h(tau_grid));

m = min(z_gamma_z);
M = max(z_gamma_z);

% figure;
% loglog(tau_grid_real, z_gamma_z)
end 